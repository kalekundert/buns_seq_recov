#!/usr/bin/env python3

import os
import warnings
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from dataclasses import dataclass
from color_me import ucsf

from .cli import main, log
from .layout import Benchmark

aa_table = {
        'ALA': 'A',
        'CYS': 'C',
        'ASP': 'D',
        'GLU': 'E',
        'GLY': 'G',
        'HIS': 'H',
        'ILE': 'I',
        'LYS': 'K',
        'LEU': 'L',
        'MET': 'M',
        'ASN': 'N',
        'PHE': 'F',
        'PRO': 'P',
        'GLN': 'Q',
        'ARG': 'R',
        'SER': 'S',
        'THR': 'T',
        'VAL': 'V',
        'TRP': 'W',
        'TYR': 'Y',
}
aa_types = {
        'D': 'positive',
        'E': 'positive',

        'K': 'negative',
        'R': 'negative',
        'H': 'negative',

        'N': 'polar',
        'Q': 'polar',
        'S': 'polar',
        'T': 'polar',

        'A': 'nonpolar',
        'V': 'nonpolar',
        'I': 'nonpolar',
        'L': 'nonpolar',
        'C': 'nonpolar',
        'M': 'nonpolar',

        'F': 'aromatic',
        'Y': 'aromatic',
        'W': 'aromatic',

        'G': 'special',
        'P': 'special',
}
type_colors = {
        'positive': ucsf.blue[0],
        'negative': ucsf.red[0],
        'polar': ucsf.purple[0],
        'nonpolar': ucsf.black[1],
        'aromatic': ucsf.black[0],
        'special': ucsf.orange[0],
}
sfxn_titles = {
        'ref': 'ref2015',
        'ref_buns_02': 'ref2015 + buried_unsatisfied_penalty (weight=0.2)',
        'ref_buns_04': 'ref2015 + buried_unsatisfied_penalty (weight=0.4)',
        'ref_buns_06': 'ref2015 + buried_unsatisfied_penalty (weight=0.6)',
        'ref_buns_08': 'ref2015 + buried_unsatisfied_penalty (weight=0.8)',
        'ref_buns_10': 'ref2015 + buried_unsatisfied_penalty (weight=1.0)',
}

# The rosetta score function is approximately in units of kcal/mol.
RT = 1.9872036e-3 * 293

def parse_score_files(jobs):
    # If the caller provided a benchmark, parse all the scores from that 
    # benchmark.
    try: jobs = jobs.all_jobs
    except AttributeError: pass

    scores = []

    for job in jobs:
        try:
            df = parse_score_file(job.outputs.score_path)
        except FileNotFoundError:
            log.warning(f"no score data found for {job.pdb} with score function {job.scorefxn}")
            continue

        df.insert(0, 'scorefxn', job.scorefxn)
        scores.append(df)

    return pd.concat(scores, ignore_index=True)

def parse_score_file(path):
    scores = []

    # Extract score data from the given file.
    with path.open() as file:
        lines = file.readlines()

    for line in lines:
        tokens = line.split()
        score = {
                'pdb': path.stem,
                'resi': int(tokens[0]),
                'wt': tokens[1],
                'mut': tokens[2],
                'score': float(tokens[3]),
        }
        scores.append(score)

    df = pd.DataFrame(scores, columns=score.keys())

    def per_residue_metrics(df):
        df['score_diff'] = df['score'] - np.min(df['score'])
        df['boltz'] = np.exp(-df['score_diff'] / RT)
        df['prob'] = df['boltz'] / np.sum(df['boltz'])
        df['log_prob'] = np.log10(df['prob'])
        df['rank'] = df['score_diff'].argsort().argsort()
        return df

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        df = df.groupby('resi').apply(per_residue_metrics)

    return df

def parse_pssms(jobs):
    counts = {}
    paths = INPUTS / 'pssm'

    for path in paths.glob('*.pssm'):
        counts.update(parse_pssm(path))

    return probs_from_counts(counts)

def parse_pssm(path):
    with path.open() as file:
        lines = file.readlines()

    pdb = path.stem
    header = lines[0].split()[2:]
    pssm = {}

    for row in lines[1:]:
        expected = {}
        resi, wt, *counts = row.split()
        resi = int(resi)
        counts = map(int, counts)

        pos = Position(
                pdb=pdb,
                resi=resi,
                wt=wt,
        )
        pssm[pos] = {}

        for mut, count in zip(header, counts):
            pssm[pos][mut] = count

    return pssm
    
def parse_hbonds(job):
    dfs = []

    for resi, wt, mut, tsv in job.outputs.hbond_paths:
        df = pd.read_csv(tsv, sep='\t')
        df.insert(0, 'pdb', job.pdb)
        df.insert(1, 'muti', resi)
        df.insert(2, 'wt', wt)
        df.insert(3, 'mut', mut)
        dfs.append(df)

    return pd.concat(dfs)


def compare_wt_recov(df, sfxn_a, sfxn_b):
    wt = df[df.mut == df.wt]
    rows = []

    for (pdb, resi), group in wt.groupby(['pdb', 'resi']):
        a = group[group.scorefxn == sfxn_a]
        b = group[group.scorefxn == sfxn_b]

        if len(a) == 0:
            log.warning(f"residue {resi} in {pdb} has scores for {sfxn_b} but not {sfxn_a}, skipping...")
            continue
        if len(b) == 0:
            log.warning(f"residue {resi} in {pdb} has scores for {sfxn_a} but not {sfxn_b}, skipping...")
            continue
        if len(a) != 1:
            raise ValueError(f"found {len(a)} scores for residue {resi} in {pdb}")
        if len(b) != 1:
            raise ValueError(f"found {len(b)} scores for residue {resi} in {pdb}")

        a = a.iloc[0]
        b = b.iloc[0]

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            row = {
                    'pdb': pdb,
                    'resi': resi,
                    'wt': a['wt'],
                    sfxn_a: a['score_diff'],
                    sfxn_b: b['score_diff'],
                    'score_diff': a['score_diff'] - b['score_diff'],
                    'prob_ratio': a['prob'] / b['prob'],
                    'log_prob_ratio': a['log_prob'] - b['log_prob'],
                    'rank_diff': a['rank'] - b['rank'],
            }

        rows.append(row)

    return pd.DataFrame(rows, columns=row.keys())


@main
def cmp_sfxn_main(args):
    """\
Compare predictions made using two different score functions.

Usage:
    {PKG_PREFIX}_cmp_sfxn <directory> [<sfxn> <ref_sfxn>] [-l | -L] [-f]

Arguments:
    <directory>
        The name of the directory comprising the benchmark.

    <sfxn> and <ref_sfxn>
        The names of the score functions to compare.  If neither is given, 
        default values will be the last and first score functions, 
        respectively, listed in the benchmark configuration.

Options:
    -f --force
        Recalculate score differences and ratios from scratch, even if cached 
        values are available.

    -l --local
        Compare scores for the most recent local run, rather than the SGE runs.

    -L --all-local
        Compare scores for all the local runs, rather than the SGE runs.  
"""
    
    bench = Benchmark(args['<directory>'])
    sfxn_a = args['<sfxn>'] or bench.scorefxns[-1]
    sfxn_b = args['<ref_sfxn>'] or bench.scorefxns[0]

    if sfxn_a not in bench.scorefxns:
        raise ValueError(f"no score function '{sfxn_a}' in benchmark '{bench}'")
    if sfxn_b not in bench.scorefxns:
        raise ValueError(f"no score function '{sfxn_b}' in benchmark '{bench}'")
    if sfxn_a == sfxn_b:
        raise ValueError(f"can't compare a score function with itself: '{sfxn_a}'")

    if args['--local']:
        jobs = [bench.last_local_job]
    elif args['--all-local']:
        jobs = bench.local_jobs
    else:
        jobs = bench.sge_jobs

    tsv_path = bench.root / f"{sfxn_a}_vs_{sfxn_b}.tsv"

    if not tsv_path.exists() or args['--force']:
        df = parse_score_files(jobs)
        cmp = compare_wt_recov(df, sfxn_a, sfxn_b)

        if len(cmp) == 0:
            print("No data that can be compared.")
            return

        cmp.to_csv(tsv_path, sep='\t')

    else:
        cmp = pd.read_csv(tsv_path, sep='\t')

    print(cmp.describe())

    try: subprocess.run(['gnumeric', tsv_path])
    except KeyboardInterrupt: print()

@main
def cmp_muts_main(args):
    """\
Compare scores for all the mutations at a single position.

Usage:
    {PKG_PREFIX}_cmp_muts <directory> [<resi>]
    {PKG_PREFIX}_cmp_muts <directory> <pdb> <resi>

Arguments:
    <directory>
        The name of the directory comprising the benchmark.

    <pdb>
        The PDB ID of the position to show.  The default is to use the 
        information from the most recent local run.

    <resi>
        The residue index to show.  The default is to use the 
        index from the most recent local run, unless that run included multiple 
        indices, in which case an index must be explicitly specified.
"""
    bench = Benchmark(args['<directory>'])
    sfxns = bench.scorefxns

    if args['<pdb>']:
        pdb = args['<pdb>']
        resi = int(args['<resi>'])

    else:
        job = bench.last_local_job
        if not job:
            raise ValueError("no local jobs; must explicitly specify <pdb> and <resi>")

        pdb = job.pdb

        if args['<resi>']:
            resi = int(args['<resi>'])
        elif len(job.resis) == 1:
            resi = job.resis[0]
        elif job.resis:
            print("Please specify one of the following residue indices: {job.resis}")
            return
        else:
            print("Please specify a residue index.")
            return

    df = parse_score_files([job])
    df = df.sort_values('rank')
    print(df)

@main
def show_muts_main(args):
    """\
Display all the mutants of a single position in pymol, with the goal of 
understanding differences in score between the mutants.

Usage:
    {PKG_PREFIX}_show_muts <directory> [<resi>]
    {PKG_PREFIX}_show_muts <directory> <sfxn> <pdb> <resi>

Arguments:
    <directory>
        The name of the directory comprising the benchmark.

    <sfxn>
        The name of the score function used to generate these structures.  The 
        default is to use the score function from the most recent local run.

    <pdb>
        The PDB ID of the position to show.  The default is to use the 
        information from the most recent local run.

    <resi>
        The residue index to show.  The default is to use the 
        index from the most recent local run, unless that run included multiple 
        indices, in which case an index must be explicitly specified.

Note that the neither the structures nor the H-bond data needed by this command 
are generated during cluster runs, because they would take too much space.  The 
typical steps one would take before running this command are:

- Run a full benchmark::

    $ ssh chef
    $ {PKG_PREFIX}_qsub <dir>

- Find mutants that would be interesting to look at more closely, e.g. 
  those where the wildtype residue is better predicted without the buried_unsat 
  term than without it::

    $ {PKG_PREFIX}_cmp_sfxn <dir>

- Re-run those mutants locally, to generate structures and H-bond data for each 
  one::

    $ {PKG_PREFIX}_local <dir> <sfxn> <resi>
    
- Open the mutants in pymol::

    $ {PKG_PREFIX}_show_muts <dir>
"""

    # Parse command-line arguments:

    bench = Benchmark(args['<directory>'])

    if args['<sfxn>'] and args['<pdb>']:
        job = Job(args['<sfxn>'], args['<pdb>'])
    else:
        job = bench.last_local_job

    if args['<resi>']:
        resi = int(args['<resi>'])
    elif len(job.resis) == 1:
        resi = job.resis[0]
    elif job.resis:
        print("Please specify one of the following residue indices: {job.resis}")
        return
    else:
        print("Please specify a residue index.")
        return

    # Load H-bond data from the benchmark:

    df = parse_hbonds(job)

    # Remove atoms that don't vary between any of the models:

    def concensus(df, col):
        return len(pd.unique(df[col])) == 1

    hbonds = pd.concat([
        group
        for key, group in df.groupby(['Residue Index', 'Atom Name'])
        if not all((
                concensus(group, 'Buried?'),
                concensus(group, 'Unsatisfied?'),
                concensus(group, 'Over-satisfied?'),
        ))
    ])

    # Create selections for all the buried/exposed unsat/oversat atoms:

    def sele_from_atom(pdb, resi, atom):
        return [f'(/{pdb}///{resi}/{atom.strip()})']

    def sele_from_hbond(hbond):
        return sele_from_atom(
                hbond['pdb'],
                hbond['Residue Index'],
                hbond['Atom Name'],
        )

    def finalize_sele(sele_list):
        return ' or '.join(sele_list) if sele_list else 'none'

    buried_unsats = []
    buried_oversats = []
    exposed_unsats = []
    exposed_oversats = []

    for i, hbond in hbonds.iterrows():
        if hbond['Buried?']:
            if hbond['Unsatisfied?']:
                buried_unsats += sele_from_hbond(hbond)
            if hbond['Over-satisfied?']:
                buried_oversats += sele_from_hbond(hbond)
        else:
            if hbond['Unsatisfied?']:
                exposed_unsats += sele_from_hbond(hbond)
            if hbond['Over-satisfied?']:
                exposed_oversats += sele_from_hbond(hbond)
        
    pymol_cmd = [
            'pymol',
                '-qx',
                job.inputs.pdb_path,
                *job.outputs.pdb_paths,
                '-d', f'set cartoon_side_chain_helper, off',
                '-d', f'set sphere_scale, 0.4',
                '-d', f'select mut, resi {resi}',
                '-d', f'select env, byres all within 8 of mut',
                '-d', f'select buried_unsats, {finalize_sele(buried_unsats)}',
                '-d', f'select exposed_unsats, {finalize_sele(exposed_unsats)}',
                '-d', f'select buried_oversats, {finalize_sele(buried_oversats)}',
                '-d', f'select exposed_oversats, {finalize_sele(exposed_oversats)}',
                '-d', f'show sticks, env',
                '-d', f'show spheres, buried_unsats',
                '-d', f'show spheres, buried_oversats',
                '-d', f'hide everything, hydro',
                '-d', f'zoom env',
                '-d', f'select none',
    ]

    subprocess.run(
            pymol_cmd, 
            stdout=open(os.devnull, 'w'),
            stderr=open(os.devnull, 'w'),
            preexec_fn=os.setpgrp,
    )

@main
def wt_recov_main(args):
    pass

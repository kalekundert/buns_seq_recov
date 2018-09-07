#!/usr/bin/env python3

import sys
import docopt
import shutil
import subprocess
from .cli import main
from .layout import Benchmark, Job
from pprint import pprint

def record_rosetta_version(bench):
    git = lambda *args: unicode(subprocess.check_output(
            ['git', '-C', str(bench.rosetta_dir)] + list(args)
    ))
    with bench.rosetta_version_path.open('w') as f:
        f.write(u'Commit: {}\n'.format(git('rev-parse', 'HEAD')))
        f.write(git('status'))

def buns_seq_recov(job, resis=None):
    if job.outputs.score_path.exists():
        raise ValueError("output file '{}' already exists; refusing to overwrite.".format(job.outputs.score_path))

    cmd = [
            job.bench.rosetta_exe('buried_unsats'),
            '-in:file:s', job.inputs.pdb_path,
            '-app:sfxn', job.scorefxn,
            '-app:out:scores', job.outputs.score_path,
            '-app:out:pdbs', job.outputs.pdb_prefix,
            '-app:out:hbonds', job.outputs.hbond_prefix,
            '-app:out:save_pdbs', job.local_run,
            '-app:out:save_hbonds', job.local_run,
            '-out:mute all',
            '-out:unmute apps',
            '-out:unmute core.init',
            '-out:unmute core.pack.pack_rotamers',
    ]
    if resis:
        cmd += [
            '-app:resis', resis,
        ]

    job.outputs.mkdirs()
    subprocess.call([str(x) for x in cmd])


@main
def qsub_main(args):
    """\
Predict sequence recovery for each test case in the given benchmark.

Usage:
    {PKG_PREFIX}_qsub <directory>

Arguments:
    <directory>
        The name of directory comprising a benchmark.  More specifically, this 
        directory must contain a file called `benchmark.json` describing the 
        parameters of the benchmark (see below for a complete reference).  The 
        directory will be populated with output files as the benchmark runs.

{CONFIG_DOCS}
"""
    bench = Benchmark(args['<directory>'])
    record_rosetta_version(bench)
    subprocess.call([
            'qsub',
            '-t', '1-{}'.format(len(bench.input_combos)),
            'buns_seq_recov_sge', str(bench.root),
    ])

def sge_main():
    bench = Benchmark(sys.argv[1])
    job = Job.from_sge_task_id(bench)
    buns_seq_recov(job)

@main
def local_main(args):
    """\
Predict sequence recovery locally for the given benchmark case(s), saving 
detailed information on the structure and H-bonding network of each mutant. 

Usage:
    {PKG_PREFIX}_local <directory> <sfxn> <pdb> [<resis>]

Arguments:
    <directory>
        The name of directory comprising a benchmark.  More specifically, this 
        directory must contain a file called `benchmark.json` describing the 
        parameters of the benchmark (see below for a complete reference).  
        Detailed output files, including PDB files and H-bond tables for each 
        mutant, will be saved in the `local` subdirectory of this directory.

    <scorefxn>
        The score function to use for the simulation.  This must be one of the 
        options specified in the configuration file.

    <pdb>
        The specific 4-letter PDB code to re-run.  This must correspond to one 
        of the PDB files referenced in the configuration file.

    <resis>
        The specific residue numbers to re-run.  If not specified, all residues 
        associated with the benchmark (e.g. all interface residues) will be 
        run.  Note that you can specify residues that weren't part of the 
        original benchmark, although I don't know why you would.

{CONFIG_DOCS}
"""
    bench = Benchmark(args['<directory>'])
    job = Job(bench, args['<pdb>'], args['<sfxn>'])

    # Record the job before we actually do anything, in case something goes 
    # wrong with the JSON.  Better to crash before doing a bunch of work.
    bench.record_local_job(job)

    buns_seq_recov(job, args['<resis>'])


@main
def clear_main(args):
    """\
Delete all results from the given benchmark.

Usage:
    {PKG_PREFIX}_config <directory>
"""
    bench = Benchmark(args['<directory>'])
    for p in bench.root.glob('*'):
        if p != bench.config_path:
            if p.is_dir():
                shutil.rmtree(str(p))
            else:
                p.unlink()
    
@main
def config_main():
    """\
Display the configuration values for the indicated benchmark.

Usage:
    {PKG_PREFIX}_config <directory>

Arguments:
    <directory>
        The name of directory comprising a benchmark.  More specifically, this 
        directory must contain a file called `benchmark.json` describing the 
        parameters of the benchmark (see below for a complete reference).  The 
        final configuration values (after applying defaults) will be displayed.

{CONFIG_DOCS}
"""
    bench = Benchmark(args['<directory>'])
    print("'pdb_dir':       '{}' ({} PDBs)".format(bench.config['pdb_dir'], len(bench.pdbs)))
    print("'scorefxns':     {} ({} sfxns)".format(bench.config['scorefxns'], len(bench.scorefxns)))
    print("'rosetta_dir':   '{}'".format(bench.config['rosetta_dir']))
    print("'rosetta_build': '{}'".format(bench.config['rosetta_build']))
    

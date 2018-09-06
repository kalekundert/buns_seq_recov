#!/usr/bin/env python2

CONFIG_DOCS = """\
Configuration:
    The following options can be specified in the `benchmark.json` file.  The 
    file should contain a single dictionary, with keys corresponding to the 
    options listed below.  Defaults will apply for any options that aren't 
    specified.  However, it is required for `benchmark.json` to exist and to 
    contain a dictionary, so if you'd like to accept all the defaults, you'll 
    still need to create a file with an empty dictionary (e.g. "{}").

    "pdb_dir"  [type: str; default: "../park2016/ref"]
        The path, relative to the benchmark directory, to the directory 
        containing the PDB file to test.  This value can contain 

    "scorefxns"  [type: list; default: ["ref", "ref_buns_10"]]
        The score functions to compare in the benchmark.  The options are "ref" 
        (i.e. ref2015) or "ref_buns_XX" (e.g. ref2015 with a penalty for buried 
        unsats), where XX indicates that the `buried_unsatisfied_penalty` has a 
        weight of X.X, and can be either "02", "04", "06", "08", or "10".

    "rosetta_dir"  [type: str; default: $ROSETTA environment variable]
        The path to the root of the Rosetta installation to use for this 
        benchmark.  The given installation must include a 'source/bin' 
        directory, which means that Rosetta needs to be compiled.

    "rosetta_build"  [type: str; default: $ROSETTA_BUILD environment variable]
        The build to use, e.g. "linuxclangrelease".  The specific builds depend 
        on what platfrom you're on.  It's best to use a release build for full 
        benchmark runs.
"""

import sys
import docopt
import shutil
import subprocess
from . import Benchmark, Job
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
            '-app:out:save_pdbs', job.run_locally,
            '-app:out:save_hbonds', job.run_locally,
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


def qsub_main():
    """\
Predict sequence recovery for each test case in the given benchmark.

Usage:
    buns_seq_recov_qsub <directory>

Arguments:
    <directory>
        The name of directory comprising a benchmark.  More specifically, this 
        directory must contain a file called `benchmark.json` describing the 
        parameters of the benchmark (see below for a complete reference).  The 
        directory will be populated with output files as the benchmark runs.

{CONFIG_DOCS}
"""
    args = docopt.docopt(qsub_main.__doc__.format(**globals()))
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

def local_main():
    """\
Predict sequence recovery locally for the given benchmark case(s), saving 
detailed information on the structure and H-bonding network of each mutant. 

Usage:
    buns_seq_recov_local <directory> <sfxn> <pdb> [<resis>]

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
    args = docopt.docopt(local_main.__doc__.format(**globals()))
    bench = Benchmark(args['<directory>'])
    job = Job(bench, args['<pdb>'], args['<sfxn>'])
    buns_seq_recov(job, args['<resis>'])

def clear_main():
    """\
Delete all results from the given benchmark.

Usage:
    buns_seq_recov_config <directory>
"""
    args = docopt.docopt(config_main.__doc__.format(**globals()))
    bench = Benchmark(args['<directory>'])
    for p in bench.root.glob('*'):
        if p != bench.config_path:
            if p.is_dir():
                shutil.rmtree(str(p))
            else:
                p.unlink()
    
def config_main():
    """\
Display the configuration values for the indicated benchmark.

Usage:
    buns_seq_recov_config <directory>

Arguments:
    <directory>
        The name of directory comprising a benchmark.  More specifically, this 
        directory must contain a file called `benchmark.json` describing the 
        parameters of the benchmark (see below for a complete reference).  The 
        final configuration values (after applying defaults) will be displayed.

{CONFIG_DOCS}
"""
    args = docopt.docopt(config_main.__doc__.format(**globals()))
    bench = Benchmark(args['<directory>'])
    print "'pdb_dir':       '{}' ({} PDBs)".format(bench.config['pdb_dir'], len(bench.pdbs))
    print "'scorefxns':     {} ({} sfxns)".format(bench.config['scorefxns'], len(bench.scorefxns))
    print "'rosetta_dir':   '{}'".format(bench.config['rosetta_dir'])
    print "'rosetta_build': '{}'".format(bench.config['rosetta_build'])
    

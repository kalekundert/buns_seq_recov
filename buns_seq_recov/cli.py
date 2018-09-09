#!/usr/bin/env python3

import docopt
import functools
import logging

PKG_PREFIX = "buns_seq_recov"
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

log = logging.getLogger(PKG_PREFIX)

def main(func):

    @functools.wraps(func)
    def decorator():
        import docopt
        args = docopt.docopt(func.__doc__.format(**globals()))
        func(args)

    return decorator



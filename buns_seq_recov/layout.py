#!/usr/bin/env python2

import os
import json
import itertools
from pathlib import Path

class Benchmark:
    default_config = {
            'pdb_dir': '../park2016/pdb/ref',
            'scorefxns': ['ref', 'ref_buns_10'],
            'rosetta_dir': os.environ.get('ROSETTA', 'rosetta'),
            'rosetta_build': os.environ.get('ROSETTA_BUILD', 'default'),
    }

    class JobInput:

        def __init__(self, bench, job):
            assert bench is job.bench
            self.bench = bench
            self.job = job

        @property
        def pdb_path(self):
            return self.bench.input_pdb_dir / '{}.pdb'.format(self.job.pdb)


    class JobOutput:

        def __init__(self, bench, job):
            assert bench is job.bench
            self.bench = bench
            self.job = job

        def mkdirs(self):
            paths = [
                    self.score_path,
                    self.pdb_prefix,
                    self.hbond_prefix,
            ]
            for p in paths:
                p.parent.mkdir(parents=True)

        @property
        def output_dir(self):
            if self.job.run_locally:
                return self.bench.root / 'local' / self.job.scorefxn 
            else:
                return self.bench.root / self.job.scorefxn 

        @property
        def score_path(self):
            return self.output_dir / 'score' / '{}.tsv'.format(self.job.pdb)

        @property
        def pdb_prefix(self):
            return self.output_dir / 'pdb' / '{}_'.format(self.job.pdb)

        @property
        def hbond_prefix(self):
            return self.output_dir / 'hbond' / '{}_'.format(self.job.pdb)

    def __init__(self, root):
        if not Path(root).is_dir():
            raise ConfigError("no such directory '{}'".format(root))

        self.root = Path(root).resolve()
        self.root_str = root  # For use in `__str__()`
        self.config_path = self.root / 'benchmark.json'
        self.config = self.default_config.copy()

        if not self.config_path.exists():
            raise ConfigError("no '{}' config file in directory '{}'".format(self.config_path.name, root))

        with self.config_path.open() as f:
            custom_config = json.load(f)

        # Update the default configuration with values from the JSON file.
        for key, value in custom_config.items():
            if key not in self.config:
                raise ConfigError("unknown configuration option '{}'".format(key))
            self.config[key] = value

        if not self.rosetta_bin_dir.exists():
            raise ConfigError("no Rosetta 'bin/' directory found: '{}'".format(self.rosetta_bin_dir))
    
    def __str__(self):
        return self.root_str

    @property
    def pdbs(self):
        return [x.stem for x in self.input_pdb_dir.glob('*.pdb')]

    @property
    def scorefxns(self):
        return self.config['scorefxns']

    @property
    def input_combos(self):
        return list(itertools.product(self.pdbs, self.scorefxns))

    @property
    def input_pdb_dir(self):
        return self.root / Path(self.config['pdb_dir'])

    def inputs(self, job):
        return self.JobInput(self, job)

    def outputs(self, job):
        return self.JobOutput(self, job)

    def rosetta_exe(self, exe):
        return self.rosetta_bin_dir / '{}.{}'.format(exe, self.rosetta_build)

    @property
    def rosetta_dir(self):
        return Path(self.config['rosetta_dir'])

    @property
    def rosetta_bin_dir(self):
        return self.rosetta_dir / 'source' / 'bin'

    @property
    def rosetta_build(self):
        return self.config['rosetta_build']
    @property
    def rosetta_version_path(self):
        return self.root / 'rosetta_version'

    @property
    def record_local_job(self, job):
        local_jobs = self.local_jobs
        local_jobs.append(job)
        local_jobs_json = [x.to_json() for x in local_jobs]

        with open(self.local_job_path, 'w') as f:
            f.write(local_jobs_json)

    @property
    def last_local_job(self):
        try:
            return self.local_jobs[-1]
        except IndexError:
            return None

    @property
    def local_jobs(self):
        if not self.local_job_path.exists():
            return []
        with open(self.local_job_path) as f:
            return [Job.from_json(self, x) for x in json.load(f)]

    @property
    def local_job_path(self):
        return self.root / 'last_local_job.json'


class Job:

    @classmethod
    def from_sge_task_id(cls, bench):
        id = int(os.environ['SGE_TASK_ID'])
        pdb, scorefxn = bench.input_combos[id]
        return cls(bench, pdb, scorefxn, run_locally=False)

    def to_json(self):
        return {'pdb': self.pdb, 'sfxn': self.scorefxn}

    @classmethod
    def from_json(cls, bench, json_dict):
        return cls(bench, json_dict['pdb'], json_dict['scorefxn'])

    def __init__(self, bench, pdb, scorefxn, run_locally=True):
        if pdb not in bench.pdbs:
            raise ValueError("no pdb '{}' in benchmark '{}'".format(pdb, bench))
        if scorefxn not in bench.scorefxns:
            raise ValueError("no score function '{}' in benchmark '{}'".format(scorefxn, bench))

        self.bench = bench
        self.pdb = pdb
        self.scorefxn = scorefxn
        self.inputs = bench.inputs(self)
        self.outputs = bench.outputs(self)
        self.run_locally = run_locally


class ConfigError(Exception):
    pass

#!/usr/bin/env python3

import os
import json
import shutil
import itertools
from pathlib import Path

class Benchmark:
    default_config = {
            'pdb_dir': '../park2016/pdb/ref',
            'scorefxns': ['ref', 'ref_buns_10'],
            'algorithm': ['pack'],
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
                p.parent.mkdir(parents=True, exist_ok=True)

        @property
        def output_dir(self):
            if self.job.local_run:
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
    
    def __eq__(self, other):
        return self.root == other.root
    
    def __str__(self):
        return self.root_str

    def clear(self):
        for p in self.root.glob('*'):
            if p != self.config_path:
                if p.is_dir():
                    shutil.rmtree(str(p))
                else:
                    p.unlink()

    @property
    def pdbs(self):
        return [x.stem for x in self.input_pdb_dir.glob('*.pdb')]

    @property
    def scorefxns(self):
        return self.config['scorefxns']

    @property
    def algorithm(self):
        return self.config['algorithm']

    def inputs(self, job):
        return self.JobInput(self, job)

    @property
    def input_combos(self):
        return list(itertools.product(self.pdbs, self.scorefxns))

    @property
    def input_pdb_dir(self):
        pdb_dir = Path(self.config['pdb_dir']).expanduser()
        if pdb_dir.is_absolute():
            return pdb_dir
        else:
            return self.root / pdb_dir

    def outputs(self, job):
        return self.JobOutput(self, job)

    @property
    def log_dir(self):
        return self.root / 'log'

    def rosetta_exe(self, exe):
        return self.rosetta_bin_dir / '{}.{}'.format(exe, self.rosetta_build)

    @property
    def rosetta_dir(self):
        rosetta_dir = Path(self.config['rosetta_dir']).expanduser()
        if rosetta_dir.is_absolute():
            return rosetta_dir
        else:
            return self.root / rosetta_dir

    @property
    def rosetta_bin_dir(self):
        return self.rosetta_dir / 'source' / 'bin'

    @property
    def rosetta_build(self):
        return self.config['rosetta_build']
    @property
    def rosetta_version_path(self):
        return self.root / 'rosetta_version'

    def jobs(self, local_run=True):
        return self.local_jobs if local_run else self.sge_jobs

    @property
    def sge_jobs(self):
        return [
                Job(self, pdb, sfxn, local_run=False)
                for pdb, sfxn in self.input_combos
        ]

    @property
    def local_jobs(self):
        if not self.local_job_path.exists():
            return []
        with open(self.local_job_path) as f:
            return [Job.from_json(self, x) for x in json.load(f)]

    @property
    def local_job_path(self):
        return self.root / 'local_jobs.json'

    @property
    def last_local_job(self):
        try:
            return self.local_jobs[-1]
        except IndexError:
            return None

    def record_local_job(self, job):
        local_jobs = self.local_jobs
        if job in local_jobs:
            return

        local_jobs.append(job)
        local_jobs_json = [x.to_json() for x in local_jobs]

        with open(self.local_job_path, 'w') as f:
            json.dump(local_jobs_json, f)


class Job:

    @classmethod
    def from_sge_task_id(cls, bench):
        id = int(os.environ['SGE_TASK_ID'])
        pdb, scorefxn = bench.input_combos[id - 1]
        return cls(bench, pdb, scorefxn, local_run=False)

    def to_json(self):
        return {'pdb': self.pdb, 'sfxn': self.scorefxn, 'resis': self.resis}

    @classmethod
    def from_json(cls, bench, json_dict):
        return cls(bench, json_dict['pdb'], json_dict['sfxn'], resis=json_dict['resis'])

    def __init__(self, bench, pdb, scorefxn, *, resis=None, local_run=True):
        if pdb not in bench.pdbs:
            raise ValueError("no pdb '{}' in benchmark '{}'".format(pdb, bench))
        if scorefxn not in bench.scorefxns:
            raise ValueError("no score function '{}' in benchmark '{}'".format(scorefxn, bench))

        def resis_from_str_list_none(resis):
            if isinstance(resis, str):
                return [int(x) for x in resis.split(',')]
            else:
                return resis or []

        self.bench = bench
        self.pdb = pdb
        self.scorefxn = scorefxn
        self.resis = resis_from_str_list_none(resis)
        self.inputs = bench.inputs(self)
        self.outputs = bench.outputs(self)
        self.local_run = local_run

    def __eq__(self, other):
        as_tuple = lambda x: (x.bench, x.pdb, x.scorefxn, x.resis, x.local_run)
        return as_tuple(self) == as_tuple(other)

    def __str__(self):
        return f"Job({self.bench}, {self.pdb}, {self.scorefxn})"

    @property
    def resis_str(self):
        return ','.join(map(str, self.resis))


class ConfigError(Exception):
    pass


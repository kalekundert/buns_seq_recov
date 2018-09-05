#!/usr/bin/env python3
# encoding: utf-8

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import re
with open('buns_seq_recov/__init__.py') as file:
    version_pattern = re.compile("__version__ = '(.*)'")
    version = version_pattern.search(file.read()).group(1)
with open('README.rst') as file:
    readme = file.read()

setup(
    name='buns_seq_recov',
    version=version,
    author='Kale Kundert',
    author_email='kale@thekunderts.net',
    long_description=readme,
    packages=[
        'buns_seq_recov',
    ],
    entry_points={
        'console_scripts': [
            'buns_seq_recov=buns_seq_recov:main',
        ],
    },
    install_requires=[
    ],
)

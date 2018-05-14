#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='kinda',
    version='0.1',
    description='Kinetic DNA strand-displacement Analyzer',
    long_description=readme,
    url='',
    author='Joseph Berleant, Chris Berlind, Erik Winfree',
    author_email='winfree@caltech.edu',
    license=license,
    classifiers=['Programming Language :: Python :: 2'],
    test_suite='tests',
    install_requires=[
      'peppercornenumerator==0.6',
      'numpy'],
    dependency_links=['http://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator/tarball/master#egg=peppercornenumerator-0.6'],
    packages=find_packages(),
    scripts=['scripts/KinDA']
)


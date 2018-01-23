# KinDA: Kinetic DNA strand-displacement Analyzer 

This package provides a framework for analyzing the kinetic behavior
of domain-level strand displacement (DSD) reaction networks with
nucleotide sequences assigned to each domain. KinDA optionally performs
domain-level reaction enumeration using the Peppercorn enumerator
(https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator).
Reaction enumeration may be skipped if detailed and condensed reactions are
given directly to KinDA. Thermodynamic and kinetic statistics about the
behavior of resting macrostates and condensed reactions are collected using
Multistrand (https://github.com/DNA-and-Natural-Algorithms-Group/multistrand),
and may be computed at a desired level of precision using KinDA.

## Dependencies
Prior to installing KinDA, make sure you have the following packages installed:
* Multistrand (https://github.com/DNA-and-Natural-Algorithms-Group/multistrand)
* Nupack

In addition, ensure the environment variable `NUPACKHOME` is set to the base directory
of your Nupack files.

KinDA will automatically install the following packages, if necessary:
* Peppercorn enumerator (https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator)


## Installation
```bash
$ python setup.py install
```
or
```
$ python setup.py install --user
```

## Quickstart using the example script "analyze.py"

### Quickstart
Load the file `Zhang_etal_Science2007.pil` in an interactive KinDA session, and 
begin querying basic data.

```sh
$ python -i analyze.py Zhang_etal_Science2007.pil
```

### Input format

PIL Files may be input in either old-style or new-style (kernel) notation.
This example file is specified in old-style notation:

```
# Zhang_etal_Science2007.pil
#
# This PIL file represents the entropy-driven catalyst system described by:
#  David Zhang, Andrew Turberfield, Bernard Yurke, Erik Winfree (Science, 2007)
#  Engineering Entropy-Driven Reactions and Networks Catalyzed by DNA
# 
# See Figure 1B for description of each system component.

# Specify all domains
sequence d1 = CTTTCCTACA : 10
sequence d2 = CCTACGTCTCCAACTAACTTACGG : 24
sequence t3 = CCCT : 4
sequence d4 = CATTCAATACCCTACG : 16
sequence t5 = TCTCCA : 6
sequence d6 = CCACATACATCATATT : 16

# Specify all strands
strand F = d2 t3 d4 : 44
strand C = d4 t5 : 22
strand OB = d1 d2 : 34
strand SB = d6 t3 d4 : 36
strand LB = t5* d4* t3* d2* : 50

# Specify all predicted complexes
structure Fuel = F : ...
structure Catalyst = C : ..
structure Substrate = OB + SB + LB : .(+.((+.)))
structure Waste = F + LB : (((+.)))
structure Signal = SB : ...
structure Output = OB : ..
structure Intermediate = OB + C + LB : .(+((+)).)
```

## Configuration Options
The file `kinda/options.py` contains optional arguments that may be modified to change the behavior of Multistrand, Peppercorn, Nupack, and KinDA.

## Version
0.5

## Authors
Joseph Berleant, Chris Berlind, and Erik Winfree



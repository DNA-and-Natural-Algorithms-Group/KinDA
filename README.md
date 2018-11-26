# KinDA: Kinetic DNA strand displacement analyzer

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
* Nupack 3.2.2+

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

### Quickstart: PIL files
Load the file `Zhang_etal_Science2007.pil` in an interactive KinDA session, and 
begin querying basic data.

For example, run the following from within the `examples` subdirectory of KinDA.

```sh
$ python -i analyze.py Zhang_etal_Science2007.pil
```

### Quickstart: Pure Python scripts
KinDA objects can be created directly in a Python script using the `kinda.objects` package. For example:
```
# simple.py
#
# This Python file creates a simple toehold-exchange system.
# The objects are identical to those described in KinDA/examples/simple.pil

# Import kinda and kinda.objects
import kinda
import kinda.objects

# Create domains
t1 = kinda.objects.Domain(name = 't1', sequence = 'AAAGAT')
d2 = kinda.objects.Domain(name = 'd2', sequence = 'AGCTGACTTA')
t3 = kinda.objects.Domain(name = 't3', sequence = 'TCCCTT')

# Create strands
strand_top1 = kinda.objects.Strand(name = 'strand_top1', domains = [t1, d2])
strand_top2 = kinda.objects.Strand(name = 'strand_top2', domains = [d2, t3])
strand_base = kinda.objects.Strand(name = 'strand_base', domains = [t3.complement, d2.complement, t1.complement])

# Create complexes
T1Bound = kinda.objects.Complex(name = 'T1Bound', strands = [strand_top1, strand_base], structure = '((+.))')
T3Intruder = kinda.objects.Complex(name = 'T3Intruder', strands = [strand_top2], structure = '..')
T3Bound = kinda.objects.Complex(name = 'T3Bound', strands = [strand_top2, strand_base], structure = '((+)).')
T1Intruder = kinda.objects.Complex(name = 'T1Intruder', strands = [strand_top1], structure = '..')

# Create System object
system = kinda.System(complexes = [T1Bound, T3Intruder, T3Bound, T1Intruder])

# Get statistics about productive condensed reactions
reactions = system.get_reactions(spurious=False, unproductive=False)
reaction = reactions[0]
system.get_stats(reaction).get_k1(relative_error = 0.5, max_sims = 500)
system.get_stats(reaction).get_k2(relative_error = 0.5, max_sims = 500)

# Similar functions can be used to get information about resting sets (i.e. resting macrostates)
```

### Input format

Currently, PIL Files must be input in old-style PIL notation. Support for new-style (kernel) notation
is planned for a future release.
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
The file `kinda/options.py` contains optional arguments that may be modified to change the default behavior of Multistrand, Peppercorn, Nupack, and KinDA. This file must be modified prior to installation to set the default behavior.

`kinda/options.py` contains four `dict` objects:
* `kinda_params`: General parameters for KinDA and its interactions with Multistrand, NUPACK, and Peppercorn
* `multistrand_params`: Parameters for Multistrand, given directly to Multistrand when constructing a Multistrand.Options object
* `nupack_params`: Parameters for NUPACK, given directly to NUPACK when calling its `sample` executable.
* `peppercorn_params`: Parameters for Peppercorn enumeration, given to the enumerator just prior to enumeration

The initialization function for the `System` object accepts an optional keyword argument for each of these `dict`s. Each keyword argument should be supplied as a `dict` object, whose key-value pairs will override the defaults in `kinda/options.py`. 

## Version
0.2

## Authors
Joseph Berleant, Chris Berlind, Stefan Badelt, Frits Dannenberg, Joseph Schaeffer, and Erik Winfree



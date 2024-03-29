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

The principles underlying KinDA are introduced in the paper
"Automated sequence-level analysis of kinetics and thermodynamics for domain-level DNA strand-displacement systems",
by Joseph Berleant, Christopher Berlind, Stefan Badelt, Frits Dannenberg, Joseph Schaeffer and Erik Winfree,
Journal of The Royal Society Interface, 2018
(https://royalsocietypublishing.org/doi/full/10.1098/rsif.2018.0107).

Questions and comments should be addressed to Erik Winfree <winfree@caltech.edu>.

## Trying out KinDA (Public AWS Image)
The easiest way to test out KinDA is through the publicly available Amazon Web Services (AWS) Amazon Machine Image (AMI). This image is available to all AWS users, and can be found in the "Community AMIs" section when creating a new EC2 instance, using the search query "KinDA v0.2".  The scripts should run on a "t2.micro" instance, but we often use "c5.9xlarge" instances for serious simulations.  `matplotlib` is installed with the `Agg` backend default, so it can output files (PDF, etc) but not produce graphics interactively.

A note about AWS regions: The KinDA AMI is currently only available in AWS's four U.S. subdivisions. If you are having trouble finding this AMI, your AWS region may be set outside the U.S. Please contact the project team if you cannot switch your account's region setting and would like us to copy the image to a new region.

## Setup and Installation
The following instructions are intended for users wishing to setup KinDA on their own machine. They should not be necessary if using the public AWS AMI, which has KinDA and all its dependencies pre-installed.

### Dependencies
KinDA requires Python 3.9+ to run.

Prior to installing KinDA, make sure you have the following packages installed:
* NUPACK 4.0.1+ (http://www.nupack.org)
* Multistrand 2.2+
  (https://github.com/DNA-and-Natural-Algorithms-Group/multistrand)

KinDA will automatically install the following packages, if necessary:
* Peppercorn enumerator 1.1+
  (https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator)

To run the plotting examples in the case studies, you will need `matplotlib`.

### Installation
```bash
$ pip install .
```
or
```
$ pip install -e .
```

## Getting Started

### Quickstart: PIL files
PIL files describe the sequences, strands, and complexes in a DNA
strand-displacement system. The file `examples/Zhang_etal_Science2007.pil`
describes an entropy-driven catalytic cascade.

The script `analyze.py` found in the `examples` subdirectory of KinDA shows how
to query basic data from a system described by a PIL file. For example, run the
following from within the `examples` subdirectory of KinDA. This will take a few
hours on a `t2.micro` instance, and proportionally less time on a faster
multiprocessor instance.

```sh
$ python -i analyze.py Zhang_etal_Science2007.pil
```

If you just want to see a script run, but don't have much time, try the (still
not so fast) simple toehold-mediated strand displacement example below, which
should take less than 100 sec of wall-clock time on a machine with 16 cores.

```sh
$ python -i analyze.py simple.pil
```

You can make all this quicker (or slower) by changing the accuracy target
(relative error) and sampling budget (number of Nupack samples / Multistrand
trajectories) at the following lines in `analyze.py`. The comments are hopefully
self-explanatory.

```
rel_error = 0.25
max_sims = 1000
```

Either way, these scripts will dump you in the python shell at the end, where
you can examine your data further, or look at the KinDA documentation, e.g.,

```
help(rxn_stats)
```

### Gentle introduction: Pure Python scripts
KinDA objects can be created directly in a Python script using the
`kinda.objects` package. For example:

```
# simple.py
#
# This Python file creates a simple toehold-exchange system.
# The objects are identical to those described in KinDA/examples/simple.pil

# Import kinda and kinda.objects
import kinda
import kinda.objects

# Create domains
t1 = kinda.objects.Domain(name='t1', sequence='AAAGAT')
d2 = kinda.objects.Domain(name='d2', sequence='AGCTGACTTA')
t3 = kinda.objects.Domain(name='t3', sequence='TCCCTT')

# Create strands
strand_top1 = kinda.objects.Strand(name='strand_top1', domains=[t1, d2])
strand_top2 = kinda.objects.Strand(name='strand_top2', domains=[d2, t3])
strand_base = kinda.objects.Strand(
    name='strand_base', domains=[t3.complement, d2.complement, t1.complement])

# Create complexes
T1Bound = kinda.objects.Complex(
    name='T1Bound', strands=[strand_top1, strand_base], structure='((+.))')
T3Intruder = kinda.objects.Complex(
    name='T3Intruder', strands=[strand_top2], structure='..')
T3Bound = kinda.objects.Complex(
    name='T3Bound', strands=[strand_top2, strand_base], structure='((+)).')
T1Intruder = kinda.objects.Complex(
    name='T1Intruder', strands=[strand_top1], structure='..')

# Create System object
system = kinda.System(complexes=[T1Bound, T3Intruder, T3Bound, T1Intruder])

# Get statistics about productive condensed reactions
print()
reactions = system.get_reactions(spurious=False, unproductive=False)
reaction = reactions[0]
system.get_stats(reaction).get_k1(relative_error=0.25, max_sims=2000, verbose=1)
system.get_stats(reaction).get_k2(relative_error=0.25, max_sims=2000, verbose=1)

# Similar functions can be used to get information about resting sets (i.e.
# resting macrostates)
```

### Gentle introduction: Input format
Currently, PIL Files must be input in old-style PIL notation. Support for
new-style (kernel) notation is planned for a future release. This example file
is specified in old-style notation:

```
# examples/Zhang_etal_Science2007.pil
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

### A deeper dive: The case studies
The directory `case_studies` contains scripts used to run the simulations
described in the paper. You should take a look at the scripts themselves, as
instructions for how to configure and run them are often included near the top.

Note that the simulations for Figure 9 were performed using a commandline
interface to KinDA (placed in the `src/kinda/scripts` directory) that is
explained in the Figure 9 directory's `README.md`.

Also note that some scripts produce plots using matplotlib, which on AWS should
be able to produce PDF output files, but you won't be able to look at graphics
interactively.

In the Amazon AMI, but not the GitHub repository, each case study directory has
a subdirectory `publication_data` with KinDA data files that were generated for
the paper. If you are not using the AMI, this data can be obtained from
(http://www.dna.caltech.edu/SupplementaryMaterial/KinDA_paper_data/).

## Configuration Options
The file `src/kinda/options.py` contains optional arguments that may be modified
to change the default behavior of Multistrand, Peppercorn, NUPACK, and KinDA.
This file must be modified prior to installation to set the default behavior,
unless the installation is in editable mode (via `pip install -e`).

`src/kinda/options.py` contains four `dict` objects:
* `kinda_params`: General parameters for KinDA and its interactions with
  Multistrand, NUPACK, and Peppercorn.
* `multistrand_params`: Parameters for Multistrand, used when constructing a
  `multistrand.Options()` object.
* `nupack_params`: Parameters for NUPACK, given directly to NUPACK when calling
  its `sample` executable.
* `peppercorn_params`: Parameters for Peppercorn enumeration, given to the
  enumerator just prior to enumeration.

The initialization function for the `kinda.System()` object accepts an optional
keyword argument for each of these `dict`s. Each keyword argument should be
supplied as a `dict` object, whose key-value pairs will override the defaults in
`src/kinda/options.py`.

## Version History
### 0.3 (August 2023)
- Migrated to Python 3.9+ and updated the Python package definition.
- Updated dependencies: Multistrand 2.2, NUPACK 4.0.1, Peppercorn enumerator
  1.1, DSDobjects 0.8.
- Created an [Apptainer](https://apptainer.org/) container for a fully
  reproducible installation.
- Reworked the order semantics of the internal object hierarchy, leading to more
  consistent analysis outputs.
- Added an optional `rate_model` key to the `multistrand_params` configuration,
  which invokes a kinetic parameter preset in Multistrand.
- Improved the reliability of `simulation.*job` by switching from the standard
  library module `multiprocessing` to the `multiprocess` package.
- Improved argument handling and output formatting in analysis scripts.

### 0.2 (April 2019)

### 0.1 (February 2018)

## Contributors
### Original authors
Joseph Berleant, Chris Berlind, Stefan Badelt, Frits Dannenberg, Joseph
Schaeffer, and Erik Winfree

### Later contributors
Boyan Beronov

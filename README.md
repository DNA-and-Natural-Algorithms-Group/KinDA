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

The principles underlying KinDA are introduced in the paper:
```
"Automated sequence-level analysis of kinetics and thermodynamics for domain-level DNA strand-displacement systems",
Joseph Berleant, Christopher Berlind, Stefan Badelt, Frits Dannenberg, Joseph Schaeffer and Erik Winfree.
Journal of The Royal Society Interface, 2018
(https://royalsocietypublishing.org/doi/full/10.1098/rsif.2018.0107).
```

General questions and comments should be addressed to Erik Winfree <winfree@caltech.edu>. For software issues regarding the most recent version,
please contact the current maintainer Boyan Beronov <beronov@cs.ubc.ca>.

## Trying out KinDA
### Public AWS Image
The easiest way to test out KinDA is through the publicly available Amazon Web Services (AWS) Amazon Machine Image (AMI). This image is available to all AWS users, and can be found in the "Community AMIs" section when creating a new EC2 instance, using the search query "KinDA v0.2".  The scripts should run on a "t2.micro" instance, but we often use "c5.9xlarge" instances for serious simulations.  `matplotlib` is installed with the `Agg` backend default, so it can output files (PDF, etc) but not produce graphics interactively.

A note about AWS regions: The KinDA AMI is currently only available in AWS's four U.S. subdivisions. If you are having trouble finding this AMI, your AWS region may be set outside the U.S. Please contact the project team if you cannot switch your account's region setting and would like us to copy the image to a new region.

### Local installation
The following instructions are intended for users wishing to setup KinDA on their
own machine. This should not be necessary if using the public AWS AMI, which has
KinDA 0.2 and all of its dependencies pre-installed.

KinDA 0.3+ runs on Python 3.10+, and requires a manual installation of the
following packages:
* NUPACK 4.0.1+ (http://www.nupack.org)
* Multistrand 2.2+
  (https://github.com/DNA-and-Natural-Algorithms-Group/multistrand)

Given the above, you can simply execute
```bash
$ pip install .
```
or
```bash
$ pip install -e .
```
in the root directoy. This will also install all remaining dependencies
from the Python Package Index, including:
* Peppercorn enumerator 1.1+
  (https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator)

### [Apptainer](https://apptainer.org/) container
Alternatively, a fully reproducible, isolated and relocatable installation
can be achieved using containerization.

1. Follow the instructions for building a Multistrand container
   (https://github.com/DNA-and-Natural-Algorithms-Group/multistrand?tab=readme-ov-file#apptainer-container).
2. Build a KinDA container layer on top, see `scripts/kinda.def`.

## Getting Started

### Gentle introduction
#### Input format
In order to use KinDA for analyzing reaction statistics,
one has to first create a `kinda.System` object, which has several functions:

1. it describes the sequences, strands and complexes in a DNA
   strand-displacement system,
2. it provides convenient access to the reactions and resting sets of such
   a system, and
3. it holds the corresponding `Stats` objects.

The user needs to provide (1), while (2) and (3) are computed by KinDA.
In particular, there are several ways of creating a reaction system description:

- Directly build up the Python objects defined in the `kinda.objects`
  subpackage, as demonstrated in `example/analyze.py::create_system()`.
- Import from an old-style PIL file, as shown in
  `example/analyze.py::import_system()`.
- Import from related Python packages such as Peppercorn or Multistrand,
  using the corresponding modules `kinda.objects.io_*`.

#### Analysis
The script `examples/analyze.py` shows how to query basic data for a system
described by a PIL file. For a comprehensive example, run the
following from within the `examples` directory, which will take a few
hours on an AWS `t2.micro` instance, and proportionally less time on a faster
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
trajectories) at the following lines in `analyze.py`. The comments are
hopefully self-explanatory.

```
rel_error = 0.25
max_sims = 1000
```

Either way, running these scripts with the `python -i` flag will dump you into
the Python shell at the end, where you can examine your data further, or look
at the KinDA documentation, e.g.,

```
help(rxn_stats)
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

### Configuration
The file `src/kinda/options.py` contains optional arguments that may be
modified to change the *default* behavior of Multistrand, Peppercorn, NUPACK,
and KinDA. This file must be modified prior to installation to set the default
behavior, unless the installation is in editable mode (via `pip install -e`).

`src/kinda/options.py` contains four `dict` objects:
* `kinda_params`: General parameters for KinDA and its interactions with
  Multistrand, NUPACK, and Peppercorn.
* `multistrand_params`: Parameters for Multistrand, used when constructing a
  `multistrand.Options()` object.
* `nupack_params`: Parameters for NUPACK, given directly to NUPACK when calling
  its `sample` executable.
* `peppercorn_params`: Parameters for Peppercorn enumeration, given to the
  enumerator just prior to enumeration.

The initialization function for the `kinda.System()` object accepts
an optional keyword argument for each of these `dict`s at runtime.
Each keyword argument should be supplied as a `dict` object, whose key-value
pairs will override the defaults.

## Version History
### 0.3 (August 2023)
- Migrated to Python 3.10+ and updated the Python package definition.
- Updated dependencies: Multistrand 2.2, NUPACK 4.0.1, Peppercorn enumerator
  1.1, DSDobjects 0.8.
- Created an [Apptainer](https://apptainer.org/) container definition for a
  fully reproducible installation, see `scripts/kinda.def`.
- Reworked the order semantics of the internal object hierarchy, leading to more
  consistent analysis outputs.
- Added an optional `rate_model` key to the `multistrand_params` configuration,
  which invokes a kinetic parameter preset in Multistrand instead of specifying
  parameters individually.
- Improved the reliability of `simulation.*job` by switching from the standard
  library module `multiprocessing` to the `multiprocess` package.
- Improved argument handling and output formatting in analysis scripts.

### 0.2 (April 2019)

### 0.1 (February 2018)

## Contributors
### Original authors
Joseph Berleant, Chris Berlind, Stefan Badelt, Frits Dannenberg, Joseph
Schaeffer, and Erik Winfree

### Current maintainer
Boyan Beronov

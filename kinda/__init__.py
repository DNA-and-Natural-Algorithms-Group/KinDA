#
#  __init__.py
#  KinDA: Kinetic DNA strand-displacement Analyzer
#

__version__ = "v0.1.6"

from kinda import System, from_pil, import_data, export_data
import statistics
import simulation
import enumeration


NUPACKERROR = """\
Attempt to interface with NUPACK 3.2.2 (Pierce lab, Caltech, www.nupack.org) failed."
NUPACKHOME environment variable must point to your NUPACK install directory.
The executables are expected to be in $NUPACKHOME/build."""

import os
if 'NUPACKHOME' not in os.environ or not os.path.isfile(
        os.path.join(os.environ['NUPACKHOME'], 'build', 'bin', 'sample')):
  raise ImportError(NUPACKERROR)


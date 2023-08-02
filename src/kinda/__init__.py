#
#  __init__.py
#  KinDA: Kinetic DNA strand-displacement Analyzer
#

from importlib import metadata

__version__ = f"v{metadata.version('kinda')}"


from .kinda import System, from_pil, import_data, export_data
from . import statistics
from . import simulation
from . import enumeration


NUPACKERROR = """\
Attempt to interface with NUPACK 3.2.2 (Pierce lab, Caltech, www.nupack.org) failed."
NUPACKHOME environment variable must point to your NUPACK source directory.
The install directory must be $NUPACKHOME/build."""

import os
if 'NUPACKHOME' not in os.environ or not os.path.isfile(
        os.path.join(os.environ['NUPACKHOME'], 'build', 'bin', 'sample')):
  raise ImportError(NUPACKERROR)


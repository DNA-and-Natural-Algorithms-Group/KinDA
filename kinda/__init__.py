#
#  __init__.py
#  KinDA: Kinetic DNA strand-displacement Analyzer
#

__version__ = "v0.1.9"

from kinda import System, from_pil, import_data, export_data
import statistics
import simulation
import enumeration

# Check that NUPACKHOME environment variable is set and 'sample' executable is available there.
# If not, print a warning message.
import os
if 'NUPACKHOME' not in os.environ or not os.path.isfile(os.path.join(os.environ['NUPACKHOME'], 'build', 'bin', 'sample')):
  print "WARNING: NUPACKHOME environment variable may not be set properly. If simulations fail, make sure NUPACKHOME is set to the directory containing your Nupack executables and that the 'sample' executable exists."
del os

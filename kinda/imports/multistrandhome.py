import os, sys

if 'MULTISTRANDHOME' in os.environ:
  multistrandhome = os.environ['MULTISTRANDHOME']
  if multistrandhome not in sys.path:
    sys.path.append(multistrandhome)
else:
  print "WARNING: Could not find the MULTISTRANDHOME environment variable. Future imports of multistrand modules may fail."
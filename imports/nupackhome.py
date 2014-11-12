import os, sys

if 'NUPACKHOME' in os.environ:
  nupackhome = os.environ['NUPACKHOME']
  if nupackhome not in sys.path:
    sys.path.append(nupackhome)
else:
  print "WARNING: Could not find the NUPACKHOME environment variable. Future imports of NUPACK modules may fail."
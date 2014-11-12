import os, sys

if 'DNAOBJECTSHOME' in os.environ:
  dnaobjectshome = os.environ['DNAOBJECTSHOME']
  if dnaobjectshome not in sys.path:
    sys.path.append(dnaobjectshome)
else:
  print "WARNING: Could not find the DNAOBJECTSHOME environment variable. Future imports of DNAObjects modules may fail."
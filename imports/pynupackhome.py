import os, sys

if 'PYNUPACKHOME' in os.environ:
  pynupackhome = os.environ['PYNUPACKHOME']
  if pynupackhome not in sys.path:
    sys.path.append(pynupackhome)
else:
  print "WARNING: Could not find the PYNUPACKHOME environment variable. Future imports of the NUPACK Python interface may fail."

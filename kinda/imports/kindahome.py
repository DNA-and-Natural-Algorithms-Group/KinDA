import os, sys

if 'KINDAHOME' in os.environ:
  kindahome = os.environ['KINDAHOME']
  if kindahome not in sys.path:
    sys.path.append(kindahome)
else:
  print "WARNING: Could not find the KINDAHOME environment variable. Future imports of KinDA modules may fail."

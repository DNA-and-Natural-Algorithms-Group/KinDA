import os, sys

## NOTE: This does NOT check that the proper enumerator is actually at the given path.
## Note that the older Wolfe enumerator also uses the ENUMERATORHOME environment variable.

if 'ENUMERATORHOME' in os.environ:
  peppercornhome = os.environ['ENUMERATORHOME']
  if peppercornhome not in sys.path:
    sys.path.append(peppercornhome)
else:
  print "WARNING: Could not find the ENUMERATORHOME environment variable. Future imports of Peppercorn modules may fail."
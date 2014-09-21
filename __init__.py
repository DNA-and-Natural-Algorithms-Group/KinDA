import os, sys
if os.environ['DNAOBJECTSHOME'] not in sys.path:
  sys.path.append(os.environ['DNAOBJECTSHOME'])
if os.environ['MULTISTRANDHOME'] not in sys.path:
  sys.path.append(os.environ['MULTISTRANDHOME'])
if os.environ['NUPACKHOME'] not in sys.path:
  sys.path.append(os.environ['NUPACKHOME'])
if os.environ['PYNUPACKHOME'] not in sys.path:
  sys.path.append(os.environ['PYNUPACKHOME'])

__all__ = ['simulation']

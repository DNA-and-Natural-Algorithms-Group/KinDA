"""
restingset.py

Copyright (c) 2014 Caltech. All rights reserved.
Coded by: Joseph Berleant (jberleant@dna.caltech.edu)

This module defines a dna `resting set`, described by a set of
Complex objects.

`resting set`
    A set of dna complexes considered by some metric to be
    capable of transitioning amongst themselves on a much faster
    timescale than other interactions.
    
"""

class RestingSet(object):
  """ Represents a resting set, a list of complexes that
  may transition amongst themselves quickly. """
  
  id_counter = 0
  
  def __init__(self, *args, **kargs):
    """
    Initialization:
    
    Keyword Arguments:
    name [type=str]                  -- Name of the macrostate. An automatic
                                        name is given if this is ommitted.
    complexes [type=list of Complex] -- List of constituent complexes.
    """
    # Assign id
    self.id = RestingSet.id_counter
    RestingSet.id_counter += 1
    
    self._object_type = 'resting-set'
    
    # Assign name
    if 'name' in kargs: self.name = kargs['name']
    else: self.name = 'resting_set{0}'.format(self.id)
    
    # Assign complexes
    self._complexes = frozenset(kargs['complexes'])
    
  @property
  def complexes(self):
    """ Returns the set of Complex objects that define this RestingSet. """
    return self._complexes
    
  def __contains__(self, complex):
    """ Returns True if the given complex is in this resting set. """
    return complex in self._complexes
    
  def __eq__(self, other):
    return self._complexes == other._complexes
  def __ne__(self, other):
    return not self.__eq__(other)
    
  def __hash__(self):
    return hash(self._complexes)
    
  def __str__(self):
    return self.__repr__()
  def __repr__(self):
    return "[" + ", ".join([repr(c) for c in self._complexes]) + "]"
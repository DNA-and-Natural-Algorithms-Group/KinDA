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

from functools import total_ordering

from .complex import Complex


@total_ordering
class RestingSet:
  """
  Represents a resting set, a list of complexes that may transition amongst
  themselves quickly.
  """
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
    
    # Assign name
    if 'name' in kargs: self.name = kargs['name']
    else: self.name = 'resting_set{0}'.format(self.id)
    
    # Assign complexes
    complexes = kargs['complexes']
    assert all(isinstance(c, Complex) for c in complexes)
    self._complexes = tuple(sorted(complexes))

  @property
  def complexes(self):
    """
    Returns a list of the Complex objects that define this RestingSet.
    """
    return list(self._complexes)

  @property
  def strands(self):
    return self._complexes[0].strands

  @property
  def sequence(self):
    return [s.sequence for s in self.strands]
    
  def __contains__(self, complex):
    """
    Returns True if the given complex is in this resting set.
    """
    return complex in self._complexes
    
  def __eq__(self, other):
    assert isinstance(other, self.__class__)
    return self._complexes.__eq__(other._complexes)

  def __lt__(self, other):
    assert isinstance(other, self.__class__)
    return self._complexes.__lt__(other._complexes)

  def __hash__(self):
    return hash(self._complexes)
    
  def __str__(self):
    return "[" + ", ".join(map(repr, self._complexes)) + "]"

  def __repr__(self):
    return str(self)

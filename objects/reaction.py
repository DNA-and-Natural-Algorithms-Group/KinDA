"""
reaction.py

Copyright (c) 2010-2014 Caltech. All rights reserved.
Coded by: Joseph Berleant (jberleant@dna.caltech.edu)

This module defines an object representing a `reaction`, either
from complex(es) to complex(es) or from resting set(s) to resting
set(s).

`reaction`
    A set of complex(es) called the 'reactants' and a set of complex(es) called
    the 'products', such that the reactants convert to the products under some
    arbitrary regime.
`resting-set reaction`
    A reaction between resting sets rather than complexes.
"""

class Reaction(object):
  """ A Reaction object is defined by a set of complexes (the 'reactants') and
  a set of complexes (the 'products')."""
  
  id_counter = 0
  
  def __init__(self, *args, **kargs):
    """
    Initialization:
  
    Keyword Arguments:
    name [type=str]
    reactants [type=Complex,required]
    products [type=Complex,required]
    """
    ## Assign unique id and increment global id counter
    self.id = Reaction.id_counter
    Reaction.id_counter += 1
    
    ## Assign object type
    self._object_type = 'reaction'
    
    ## Assign object name if supplied, or create an automatic one
    if 'name' in kargs: self.name = kargs['name']
    else: self.name = 'reaction_' + self.id
    
    ## Assign reactants and products
    self._reactants = frozenset(kargs['reactants'])
    self._products = frozenset(kargs['products'])
    
    ## Optional modifiers
    self.modifiers = {}
  
  # Basic properties
  @property
  def reactants(self):
    return self._reactants
  @property
  def products(self):
    return self._products
    
  # Basic queries
  def is_reactant(self, c):
    return c in self._reactants
  def is_product(self, c):
    return c in self._products
    
  # Equality definition
  def __eq__(self, other):
    return (self._reactants == other._reactants
        and self._products == other._products)
  def __ne__(self, other):
    return not self.__eq__(other)

  def __hash__(self, other):
    return hash((self._reactants, self._products))
    
    
  
class RestingSetReaction(Reaction):
  """ The RestingSetReaction class is ostensibly for only reactions between
  resting sets. However, it is currently identical to the Reaction class and
  the two are interchangeable."""
  pass
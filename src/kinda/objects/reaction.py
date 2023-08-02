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
    reactants [type=iterable of Complex,required]
    products [type=iterable of Complex,required]
    """
    ## Assign unique id and increment global id counter
    self.id = Reaction.id_counter
    Reaction.id_counter += 1
    
    ## Assign object type
    self._object_type = 'reaction'
    
    ## Assign object name if supplied, or create an automatic one
    if 'name' in kargs: self.name = kargs['name']
    else: self.name = 'reaction_' + str(self.id)
    
    ## Assign reactants and products
    self._reactants = tuple(sorted(kargs['reactants'], key = lambda c: c.id))
    self._products = tuple(sorted(kargs['products'], key = lambda c: c.id))
    
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

  def has_reactants(self, complexes):
    return all([self._reactants.count(c) >= complexes.count(c) for c in complexes])

  def reactants_equal(self, complexes):
    reactants = list(self._reactants)
    for c in complexes:
      if reactants.count(c) == 0:
        return False
      reactants.remove(c)
    return reactants == []

  def is_product(self, c):
    return c in self._products

  def has_products(self, complexes):
    # kinda inefficient...
    return all([self._products.count(c) >= complexes.count(c) for c in complexes])

  def products_equal(self, complexes):
    products = list(self._products)
    for c in complexes:
      if products.count(c) == 0:
        return False
      products.remove(c)
    return products == []
    
  # Equality definition
  def __eq__(self, other):
    return (self._reactants == other._reactants
        and self._products == other._products)
  def __ne__(self, other):
    return not self.__eq__(other)

  def __hash__(self):
    return hash((self._reactants, self._products))
    
  # Output
  def __str__(self):
    return self.__repr__()

  def __repr__(self):
    reactant_str = ' + '.join([repr(r) for r in sorted(self._reactants, key=lambda c:c.name)])
    product_str = ' + '.join([repr(p) for p in sorted(self._products, key=lambda c:c.name)])
    return reactant_str + " -> " + product_str
    
    
  
class RestingSetReaction(Reaction):
  """
  The RestingSetReaction class is ostensibly for only reactions between resting
  sets. However, it is currently identical to the Reaction class and the two are
  interchangeable.
  """
  def __repr__(self):
    reactant_str = ' + '.join([repr(r) for r in sorted(self._reactants, key=lambda rs:[c.name for c in rs.complexes])])
    product_str = ' + '.join([repr(p) for p in sorted(self._products, key=lambda rs:[c.name for c in rs.complexes])])
    return reactant_str + " -> " + product_str

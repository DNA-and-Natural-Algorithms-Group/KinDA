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

from functools import total_ordering

from .complex import Complex
from .restingset import RestingSet


@total_ordering
class Reaction:
  """
  A Reaction object is defined by a set of complexes (the 'reactants') and a set
  of complexes (the 'products').
  """
  agent_type = Complex
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
    
    ## Assign object name if supplied, or create an automatic one
    if 'name' in kargs: self.name = kargs['name']
    else: self.name = 'reaction_' + str(self.id)
    
    ## Assign reactants and products
    reactants = kargs['reactants']
    products = kargs['products']
    assert all(isinstance(a, self.agent_type) for a in reactants)
    assert all(isinstance(a, self.agent_type) for a in products)
    self._reactants = tuple(sorted(reactants))
    self._products = tuple(sorted(products))

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
    return (self._reactants, self._products).__eq__(
      (other._reactants, other._products))

  def __lt__(self, other):
    return (self._reactants, self._products).__lt__(
      (other._reactants, other._products))

  def __hash__(self):
    return hash((self._reactants, self._products))
    
  # Output
  def __str__(self):
    return self.__repr__()

  def __repr__(self):
    reactant = ' + '.join(map(repr, self._reactants))
    product = ' + '.join(map(repr, self._products))
    return f"{reactant} -> {product}"


@total_ordering
class RestingSetReaction(Reaction):
  """
  The RestingSetReaction class is ostensibly for only reactions between resting
  sets. However, it is currently identical to the Reaction class and the two are
  interchangeable.
  """
  agent_type = RestingSet

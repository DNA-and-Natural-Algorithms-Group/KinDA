## IMPORTS
import itertools as it

from KinDA.imports import peppercornhome, dnaobjectshome

import dnaobjects as dna
import enumerator.enumerator as enum


## GLOBALS
# The follow_fast_reactions function is copied shamelessly from KinD's
# utilities.py.
def follow_fast_reactions(complexes, all_fast_reactions, restingsets, inprogress, visited):
  """Returns a list containing lists of resting states. Each of the inner lists
  represents a set of products that can result from fast reactions starting 
  from the given complex.
  An important optimization would be to add memoization. [NOTE: memoization added but untested]"""

  for i, cmplx in enumerate(complexes):
  
    # If the complex is currently being processed, you've entered a cycle, so stop enumerating
    if cmplx in inprogress:
      return []

    # If the complex has already been processed, don't reprocess it!
    if cmplx in visited:
      continue
    
    # If the complex is in a resting state, the resting state is the only product
    rs = dna.utils.get_containing_set(restingsets, cmplx)
    if rs is not None:
      visited[cmplx] = [[rs]]
      continue
    
    # Mark the complex as being processed to avoid reaction cycles
    inprogress.add(cmplx)

    # Initialize list of possible products from this complex
    prods = []
    
    # Otherwise, check all fast reactions with cmplx as the reactant
    fast_reactions = filter(lambda x: x.is_reactant(cmplx), all_fast_reactions)
    for rxn in fast_reactions:
      # Find the product sets of each of rxn's products
      rxn_results = follow_fast_reactions(rxn.products, all_fast_reactions, restingsets, inprogress, visited)
      #rxn_results = follow_fast_reactions(rxn.products, all_fast_reactions, restingsets, inprogress, dict())
      
      # Add this reaction's product sets to the cumulative list
      prods.extend(rxn_results)

    # Unmark the complex as being processed
    inprogress.remove(cmplx)
    
    # Memoize the enumerated products
    visited[cmplx] = prods
    print len(visited), len(inprogress)
    if len(inprogress) % 25 == 0: print len(inprogress)
    
  # Calculate combinatorial fates
  complex_products = [visited[cmplx] for cmplx in complexes]
  all_products = [sum(p, []) for p in it.product(*complex_products)]
  return all_products


## CLASSES
class EnumerateJob(object):
  def __init__(self, *args, **kargs):
    # Pull out domains, strands, complexes
    self.domains = kargs.get('domains', [])
    self.strands = kargs.get('strands', [])
    self.complexes = kargs.get('complexes', [])
    
    # Flag indicating if enumeration has occurred yet
    self.enumerated = False
    self.condensed = False
    
  def enumerate(self):
    ## Convert to Peppercorn objects
    enum_objects = dna.io_Peppercorn.to_Peppercorn(
        domains = self.domains,
        strands = self.strands,
        complexes = self.complexes
    )
    
    ## Instantiate and run enumerator
    e = enum.Enumerator(
        [v for k, v in enum_objects['domains']],
        [v for k, v in enum_objects['strands']],
        [v for k, v in enum_objects['complexes']]
    )
    e.enumerate()
    
    ## Convert enumerated results back to DNAObjects objects
    dna_objects = dna.io_Peppercorn.from_Peppercorn(
        reactions = e.reactions,
        restingsets = e.resting_states
    )
    rxns_dict = dict(dna_objects['reactions'])
    self.enumerated_complexes = [v for k, v in dna_objects['complexes']]
    self.enumerated_slow_rxns = set(sum(
        [[rxns_dict[rxn] for rxn in e.get_slow_reactions(c) if rxn in rxns_dict]
          for c
          in e.resting_complexes],
        []))
    self.enumerated_fast_rxns = set(sum(
        [[rxns_dict[rxn] for rxn in e.get_fast_reactions(c) if rxn in rxns_dict]
          for c
          in e.complexes],
        []))
    self.enumerated_restingsets = [v for k, v in dna_objects['restingsets']]
  
    self.enumerated = True
    
  def condense_reactions(self):  
    slow_rxns = self.enumerated_slow_rxns
    fast_rxns = self.enumerated_fast_rxns
    restingsets = self.enumerated_restingsets
    
    condensed_rxns = set()
    for rxn in slow_rxns:
      products = follow_fast_reactions(
          rxn.products,
          fast_rxns,
          restingsets,
          set(),
          dict()
      )
      reactants = [dna.utils.get_containing_set(self.enumerated_restingsets, c)
            for c
            in rxn.reactants]
      
      for p in products:
        rs_rxn = dna.RestingSetReaction(reactants = reactants, products = p)
        condensed_rxns.add(rs_rxn)

    # Make sure the reverse reaction between every pair of reactants is included
    # This is an important difference between our enumeration and Peppercorn enumeration,
    # which may not include the reverse reaction
    reactant_pairs = it.product(restingsets, restingsets)
    for reactants in reactant_pairs:
      rs_rxn = dna.RestingSetReaction(reactants = reactants, products = reactants)
      condensed_rxns.add(rs_rxn)

    self.condensed_reactions = condensed_rxns
    
  def get_restingsets(self):
    if not self.enumerated: self.enumerate()
    return self.enumerated_restingsets[:]
    
  def get_reactions(self):
    if not self.enumerated:  self.enumerate()
    return list(self.enumerated_fast_rxns) + list(self.enumerated_slow_rxns)
    
  def get_restingset_reactions(self):
    if not self.enumerated:  self.enumerate()
    if not self.condensed:  self.condense_reactions()
    return list(self.condensed_reactions)
  

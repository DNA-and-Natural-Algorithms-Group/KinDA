## IMPORTS
import itertools as it

from KinDA.imports import peppercornhome, dnaobjectshome

from .. import objects as dna
import peppercornenumerator as enum

from KinDA import options


## GLOBALS
# The follow_fast_reactions function is copied shamelessly from KinD's
# utilities.py.
def follow_fast_reactions(complexes, all_fast_reactions, restingsets, inprogress, products_dict = dict()):
  """Returns a list containing lists of resting states. Each of the inner lists
  represents a set of products that can result from fast reactions starting 
  from the given complex.
  An important optimization would be to add memoization. [NOTE: memoization added but untested]"""

  products_dict_temp = {}
  for i, cmplx in enumerate(complexes):
  
    # If the complex is currently being processed, you've entered a cycle, so stop enumerating
    if cmplx in inprogress:
      return []

    # If the complex has already been processed, don't reprocess it!
    if cmplx in products_dict_temp:
      continue
    elif cmplx in products_dict:
      products_dict_temp[cmplx] = products_dict[cmplx]
      continue
    
    # If the complex is in a resting state, the resting state is the only product
    rs = dna.utils.get_containing_set(restingsets, cmplx)
    if rs is not None:
      #visited[cmplx] = [[rs]]
      products_dict_temp[cmplx] = [[rs]]
      continue
    
    # Mark the complex as being processed to avoid reaction cycles
    inprogress.add(cmplx)

    # Initialize list of possible products from this complex
    prods = []
    
    # Otherwise, check all fast reactions with cmplx as the reactant
    fast_reactions = filter(lambda x: x.is_reactant(cmplx), all_fast_reactions)
    for rxn in fast_reactions:
      # Find the product sets of each of rxn's products
      rxn_results = follow_fast_reactions(rxn.products, all_fast_reactions, restingsets, inprogress, products_dict)
      #rxn_results = follow_fast_reactions(rxn.products, all_fast_reactions, restingsets, inprogress, dict())
      
      # Add this reaction's product sets to the cumulative list
      prods.extend(rxn_results)

    # Unmark the complex as being processed
    inprogress.remove(cmplx)

    # Store products that could result from cmplx
    products_dict_temp[cmplx] = prods
    
    # Memoize the enumerated products, but only if inprogress is empty (otherwise some paths may be unchecked)
    if len(inprogress) == 0:
      products_dict[cmplx] = prods
    #if len(inprogress) % 25 == 0: print len(inprogress)
    
  # Calculate combinatorial fates
  complex_products = [products_dict_temp[cmplx] for cmplx in complexes]
  all_products = [sum(p, []) for p in it.product(*complex_products)]
  return all_products


## CLASSES
class EnumerateJob(object):
  def __init__(self, *args, **kargs):
    # Pull out domains, strands, complexes
    self._init_reactions = kargs.get('reactions', [])
    self._init_complexes = kargs.get('complexes', [])
    self._init_strands = kargs.get('strands', [])
    self._init_domains = kargs.get('domains', [])

    self._reactions = set(self._init_reactions)
    self._complexes = set(self._init_complexes
                        + [c for rxn in self._reactions for c in rxn.reactants+rxn.products])
    self._strands = set(self._init_strands
                     + [s for c in self.complexes for s in c.strands])
    self._domains = set(self._init_domains
                     + [d for s in self.strands for d in s.base_domains()])
    
    # Flag indicating if enumeration has occurred yet
    self._enumerated = False
    self._condensed = False

    # Default vals enumerated/condensed objects
    self._enumerated_complexes = []
    self._enumerated_restingsets = []
    self._enumerated_slow_rxns = []
    self._enumerated_fast_rxns = []
    self._condensed_rxns = []
      
  @property
  def reactions(self):
    return list(self._reactions)
  @property
  def complexes(self):
    return list(self._complexes)
  @property
  def strands(self):
    return list(self._strands)
  @property
  def domains(self):
    return list(self._domains)

  @property
  def enumerated_complexes(self):
    return self._enumerated_complexes
  @property
  def enumerated_restingsets(self):
    return self._enumerated_restingsets
  @property
  def enumerated_slow_rxns(self):
    return self._enumerated_slow_rxns
  @property
  def enumerated_fast_rxns(self):
    return self._enumerated_fast_rxns
  @property
  def condensed_rxns(self):
    return self._condensed_rxns

  @property
  def enumerated(self):
    return self._enumerated
  @property
  def condensed(self):
    return self._condensed
    
  def enumerate(self):
    ## Convert to Peppercorn objects
    enum_objects = dna.io_Peppercorn.to_Peppercorn(
        domains = self.domains,
#        strands = self.strands,
        complexes = self.complexes,
        reactions = self.reactions
    )
    
    ## Instantiate and run enumerator
#    e = enum.Enumerator(
#        [v for k, v in enum_objects['domains']],
#        [v for k, v in enum_objects['strands']],
#        [v for k, v in enum_objects['complexes']]
#    )
    e = enum.Enumerator(
        initial_complexes = [v for _,v in enum_objects['complexes']],
        initial_reactions = [v for _,v in enum_objects['reactions']]
    )
    # Set peppercorn options
    if '--release-cutoff-1-1' in options.peppercorn_params:
      e.RELEASE_CUTOFF_1_1 = options.peppercorn_params['--release-cutoff-1-1']
    if '--release-cutoff-1-n' in options.peppercorn_params:
      e.RELEASE_CUTOFF_1_N = options.peppercorn_params['--release-cutoff-1-n']
    # Perform enumeration
    e.enumerate()
    
    ## Convert enumerated results back to DNAObjects objects
    dna_objects = dna.io_Peppercorn.from_Peppercorn(
        complexes = e.complexes,
        reactions = e.reactions,
        restingsets = e.resting_sets,
        domain_seqs = {d.name: str(d.constraints) for d in self.domains}
    )
    rxns_dict = dict(dna_objects['reactions'])
    self._enumerated_complexes = [v for _, v in dna_objects['complexes']]
    self._enumerated_slow_rxns = set(sum(
        [[rxns_dict[rxn] for rxn in e.get_slow_reactions(c) if rxn in rxns_dict]
          for c
          in e.resting_complexes],
        []))
    self._enumerated_fast_rxns = set(sum(
        [[rxns_dict[rxn] for rxn in e.get_fast_reactions(c) if rxn in rxns_dict]
          for c
          in e.complexes],
        []))
    self._enumerated_restingsets = [v for _, v in dna_objects['restingsets']]
  
    self._enumerated = True
    
  def condense_reactions(self):  
    slow_rxns = self.enumerated_slow_rxns
    fast_rxns = self.enumerated_fast_rxns
    restingsets = self.enumerated_restingsets
    
    condensed_rxns = set()
    i = 0
    for rxn in slow_rxns:
      products = follow_fast_reactions(
          rxn.products,
          fast_rxns,
          restingsets,
          set()#,
          #dict()
      )
      reactants = [dna.utils.get_containing_set(self.enumerated_restingsets, c)
            for c
            in rxn.reactants]
      
      for p in products:
        rs_rxn = dna.RestingSetReaction(reactants = reactants, products = p)
        condensed_rxns.add(rs_rxn)

      i += 1
      print "KinDA: Performing reaction condensation... {0}%\r".format(100*i/len(slow_rxns)),

    # Make sure the reverse reaction between every pair of reactants is included
    # This is an important difference between our enumeration and Peppercorn enumeration,
    # which may not include the reverse reaction
    reactant_pairs = it.product(restingsets, restingsets)
    for reactants in reactant_pairs:
      rs_rxn = dna.RestingSetReaction(reactants = reactants, products = reactants)
      condensed_rxns.add(rs_rxn)

    print "KinDA: Performing reaction condensation... Done!"

    self._condensed_rxns = condensed_rxns

    self._condensed = True
    
  def get_complexes(self):
    if not self.enumerated: self.enumerate()
    return self.enumerated_complexes[:]

  def get_restingsets(self):
    if not self.enumerated: self.enumerate()
    return self.enumerated_restingsets[:]
    
  def get_reactions(self):
    if not self.enumerated:  self.enumerate()
    return list(self.enumerated_fast_rxns) + list(self.enumerated_slow_rxns)
    
  def get_restingset_reactions(self):
    if not self.enumerated:  self.enumerate()
    if not self.condensed:  self.condense_reactions()
    return list(self.condensed_rxns)
  

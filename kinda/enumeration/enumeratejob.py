## IMPORTS
import itertools as it

import peppercornenumerator as enum

from .. import objects as dna
from .. import options

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
    self._enumerated_slow_reactions = []
    self._enumerated_fast_reactions = []
    self._condensed_reactions = []

    # Pull out enumeration options
    self._peppercorn_params = kargs.get('peppercorn_params', {})
      
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
  def enumerated_slow_reactions(self):
    return self._enumerated_slow_reactions
  @property
  def enumerated_fast_reactions(self):
    return self._enumerated_fast_reactions
  @property
  def condensed_reactions(self):
    return self._condensed_reactions

  @property
  def enumerated(self):
    return self._enumerated
  @property
  def condensed(self):
    return self._condensed

  @property
  def peppercorn_params(self):
    return dict(self._peppercorn_params)
    
  def enumerate(self):
    ## Convert to Peppercorn objects
    enum_objects = dna.io_Peppercorn.to_Peppercorn(
        domains = self.domains,
        complexes = self.complexes,
        reactions = self.reactions
    )
    
    ## Instantiate and run enumerator
    e = enum.Enumerator(
        initial_complexes = [v for _,v in enum_objects['complexes']],
        initial_reactions = [v for _,v in enum_objects['reactions']]
    )
    # Set peppercorn options
    for k,v in self._peppercorn_params.iteritems():
      setattr(e, k, v)

    # Perform enumeration
    print "KinDA: Performing reaction enumeration with Peppercorn...",
    e.enumerate()
    print "Done!"

    # Perform reaction condensation
    print "KinDA: Performing reaction condensation with Peppercorn...",
    enumc = enum.ReactionGraph(e)
    enumc.condense()
    print "Done!"
    
    self._enumerator = e
    self._enumerator_condensor = enumc
  
    ## Convert enumerated results back to DNAObjects objects
    domain_info = {d.name: [('seq', str(d.sequence))] for d in self.domains}
    strand_info = {tuple([d.name for d in s.base_domains()]): [('name', s.name)] for s in self.strands}
    dna_objects = dna.io_Peppercorn.from_Peppercorn(
        complexes = e.complexes,
        reactions = e.reactions,
        restingsets = e.resting_sets,
        restingsetreactions = enumc.condensed_reactions,
        domain_info = domain_info,
        strand_info = strand_info
    )
    rs_rxns_dict = dict(dna_objects['restingsetreactions'])
    rxns_dict = dict(dna_objects['reactions'])
    self._enumerated_complexes = [v for _, v in dna_objects['complexes']]
    self._enumerated_slow_reactions = set(rxns_dict[rxn] for rxn in e.reactions if rxn.rtype == 'bind21')
    self._enumerated_fast_reactions = set(rxns_dict[rxn] for rxn in e.reactions if rxn.rtype != 'bind21')
    self._enumerated_restingsets = set(v for _, v in dna_objects['restingsets'])
    self._condensed_reactions = set(v for _, v in dna_objects['restingsetreactions'])

    ## Make sure the reverse reaction between every pair of reactants is included
    ## This is an important difference between our enumeration and Peppercorn enumeration,
    ## which may not include the reverse reaction
    reactant_pairs = it.product(self._enumerated_restingsets, self._enumerated_restingsets)
    for reactants in reactant_pairs:
      rs_rxn = dna.RestingSetReaction(reactants = reactants, products = reactants)
      self._condensed_reactions.add(rs_rxn)

    self._enumerated = True
    self._condensed = True
    
   
  def get_complexes(self):
    if not self.enumerated: self.enumerate()
    return list(self.enumerated_complexes)

  def get_restingsets(self):
    if not self.enumerated: self.enumerate()
    return list(self.enumerated_restingsets)
    
  def get_reactions(self):
    if not self.enumerated:  self.enumerate()
    return list(self.enumerated_fast_reactions) + list(self.enumerated_slow_reactions)
    
  def get_restingset_reactions(self):
    if not self.enumerated or not self.condensed:  self.enumerate()
#    if not self.condensed:  self.condense_reactions()
    return list(self.condensed_reactions)
  

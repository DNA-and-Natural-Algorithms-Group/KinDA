# kinda.py
# Created by Joseph Berleant, 1/12/2018
#
# Defines the KinDA class, encapsulating statistics calculations for DNA strand-displacement system properties.


## IMPORTS

from KinDA.enumeration.enumeratejob import EnumerateJob
from KinDA.statistics import stats_utils
  

## GLOBALS

# Convenience function to create KinDA object from a given PIL file
def make_kinda_from_pil(enumerate = True):
  pass ## TODO


## CLASSES
class KinDA(object):
  """ Stores and manages Stats objects for each system component for easier retrieval.
      Data can also be stored in a file and retrieved for later analysis.
      A KinDA object is instantiated with an EnumerateJob object, from which
      detailed and condensed reactions as well as resting sets and complexes are taken. """

  def __init__(self, complexes, restingsets = [], reactions = [], rs_reactions = [], enumerate = True, c_max = None):
    """ Constructs a KinDA object with the given complexes, restingsets, reactions, and condensed reactions.
    If enumerate is True (default), the Peppercorn enumerator is used to enumerate a detailed reaction network.
    KinDA performs reaction condensation to produce the condensed reaction network. The given restingsets and
    reactions are added to the enumerated network. Detailed reactions are added to the network prior to condensation.
    If enumerate is False, then no enumeration is performed and only the given resting sets and reactions are
    used for analysis.
    c_max defines a default max concentration for each resting set, used to compute the effect of spurious and
    unproductive reactions on the effective concentrations of each complex available for intended reactions. """

    # Store given data
    # We extract any resting sets and complexes from higher-order objects to maintain consistency
    self._rs_reactions = set(rs_reactions)
    self._reactions = set(reactions)
    self._restingsets = set(restingsets
        + [rs for rxn in self._rs_reactions for rs in rxn.reactants+rxn.products]
    )
    self._complexes = set(complexes
        + [c for rs in self._restingsets for c in rs.complexes]
        + [c for rxn in self._reactions for c in rxn.reactants+rxn.products]
    )

    if enumerate:
      # Create enumeration object
      self._enum_job = EnumerateJob(complexes = self.complexes)#, reactions = reactions)### TODO: IMPLEMENT THIS FUNCTIONALITY IN EnumerateJob
      # Incorporate enumerated data
      self._complexes |= self.enum_job.get_complexes()
      self._restingsets |= self.enum_job.get_restingsets()
      self._reactions |= self.enum_job.get_reactions()
      self._rs_reactions |= self.enum_job.get_restingset_reactions()
    else:
      self._enum_job = None

    # Create stats objects for reactions and resting sets
    # make_stats() will also make stats objects for potential spurious reactions and resting sets
    self._rs_to_stats, self._rxn_to_stats = stats_utils.make_stats(self.enum_job) ### TODO: fix stats_utils.make_stats so it doesn't require enum_job

    # Pull out spurious reactions/resting sets
    self._spurious_restingsets = set(self._rs_to_stats.keys()) - self._restingsets)
    self._spurious_rs_reactions = set(self._rxn_to_stats.keys()) - self._rs_reactions)

    # Set default max concentration for each resting set
    for rs_stats in self._rs_to_stats.values():
      rs_stats.c_max = c_max


  ## Basic get functions for system objects

  @property
  def complexes(self):
    """ Returns a list of all complexes (given and enumerated) predicted for the system. """
    return list(self._complexes)

  @property
  def restingsets(self):
    """ Returns a list of all resting sets (given, enumerated, and spurious) predicted for the system. """
    return list(self._restingsets + self._spurious_restingsets)

  @property
  def reactions(self):
    """ Returns a list of detailed reactions (given and enumerated) predicted for the system. """
    return list(self._reactions)

  @property
  def restingset_reactions(self):
    """ Returns a list of all resting-set (condensed) reactions (given, enumerated, and spurious) for the system. """
    return list(self._rs_reactions + self._spurious_rs_reactions)


  ## Convenience filters for specific objects

  def get_reaction(self, reactants, products):
    """ Returns a single reaction with exactly the same reactants and products as those given. """
    rxns = self.get_reactions(reactants, products)
    rxns = filter(lambda x: x.reactants_equal(reactants) and x.products_equal(products), rxns)
    if len(rxns) >= 1:
      return next(iter(rxns))
    else:
      print "Warning: Could not find reaction with the given reactants {0} and products {1}".format(reactants, products)
      return None

  def get_reactions(self, reactants = [], products = [], unproductive = None, spurious = None):
    """ Returns a list of all reactions including the given reactants and the given products.
        If specified, spurious = True will return only spurious reactions (those not enumerated by Peppercorn)
        and spurious = False will return only enumerated reactions. Otherwise, no distinction will be made.
        """
    if spurious == True:
      rxns = self.spurious_rs_reactions
    elif spurious == False:
      rxns = self.rs_reactions
    else:
      rxns = self.spurious_rs_reactions + self.rs_reactions

    if unproductive == True:
      rxns = filter(lambda x: x.has_reactants(x.products) and x.has_products(x.reactants), rxns)
    elif unproductive == False:
      rxns = filter(lambda x: not(x.has_reactants(x.products) and x.has_products(x.reactants)), rxns)

    return filter(lambda x: x.has_reactants(reactants) and x.has_products(products), rxns)

  def get_restingsets(self, complex = None, strands = [], spurious = False):
    """ Returns a list of resting sets satisfying the filter arguments.
    If complex is specified, returns a list of resting sets containing the given complex.
    If complex is not specified and strands is specified, returns a list of resting sets containing all given strands.
    If spurious is True, only spurious resting sets are returned.
    If spurious is False, only non-spurious resting sets are returned.
    If spurious is None, both spurious and non-spurious resting sets may be returned.
    """
    if spurious == True:
      rs = self.spurious_restingsets
    elif spurious == False:
      rs = self.restingsets
    else:
      rs = self.restingsets + self.spurious_restingsets

    if complex != None:
      rs = filter(lambda x: complex in x, rs)
    else:
      rs = filter(lambda x: all([s in x.strands for s in strands]), rs)
    return rs

  def get_stats(self, obj):
    """ Returns the stats object corresponding to the given system object.
    obj must be a resting-set reaction or resting set in the system. """
    if obj in self.rxn_to_stats:
      return self.rxn_to_stats[obj]
    elif obj in self.rs_to_stats:
      return self.rs_to_stats[obj]
    else:
      print "Statistics for object {0} not found.".format(obj)



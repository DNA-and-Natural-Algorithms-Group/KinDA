## IMPORTS

from KinDA.enumeration.enumeratejob import EnumerateJob
from KinDA.statistics import stats_utils
  
## GLOBALS

## CLASSES
class KinDA(object):
  """ Stores and manages Stats objects for each system component for easier retrieval.
      Data can also be stored in a file and retrieved for later analysis.
      A SystemStats object is instantiated with an EnumerateJob object, from which
      detailed and condensed reactions as well as resting sets and complexes are taken. """

  def __init__(self, complexes, c_max = None):
    ## 
    self.enum_job = EnumerateJob(complexes = complexes)

    self.complexes = self.enum_job.get_complexes()
    self.restingsets = self.enum_job.get_restingsets()
    self.reactions = self.enum_job.get_reactions()
    self.rs_reactions = self.enum_job.get_restingset_reactions()

    self.rs_to_stats, self.rxn_to_stats = stats_utils.make_stats(self.enum_job)

    self.spurious_restingsets = list(set(self.rs_to_stats.keys()) - set(self.restingsets))
    self.spurious_rs_reactions = list(set(self.rxn_to_stats.keys()) - set(self.rs_reactions))

    for rs_stats in self.rs_to_stats.values():
      rs_stats.c_max = c_max



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
    if obj in self.rxn_to_stats:
      return self.rxn_to_stats[obj]
    elif obj in self.rs_to_stats:
      return self.rs_to_stats[obj]
    else:
      print "Statistics for object {0} not found.".format(obj)



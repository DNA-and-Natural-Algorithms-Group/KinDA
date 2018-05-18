# stats.py
# Created by Joseph Berleant, 11/11/2014
#
# Defines RestingSetRxnStats and RestingSetStats, convenience classes that encapsulate
# functions for querying statistics about DSD reactions and resting macrostates.

#### Possible future TODOs:
####   -Implement ComplexStats and ComplexRxnStats
####   -RestingSetRxnStats can be updated to take a DNAObjects.Reaction object,
####    rather than reactants and products separately.

## IMPORTS

import math

from ..simulation.multistrandjob import FirstPassageTimeModeJob, \
    TransitionModeJob, FirstStepModeJob
from ..simulation.nupackjob import NupackSampleJob
from ..enumeration.enumeratejob import EnumerateJob
from ..objects import utils, Macrostate
from .. import nupack, options
import stats_utils
    

## CLASSES

class RestingSetRxnStats(object):
  """ Calculates statistics for this resting-set reaction.
      The following statistics are calculated directly:
        - rate (via Multistrand) of collision leading to success (k1)
        - minimum k1 rate based on worst-case reactant depletion
        - rate (via Multistrand) of successful folding after collision (k2)
      The following information should be stored with this object for
      convenience in displaying and calculating certain information:
        - RestingSetStats objects for each reactant
        - ComplexRxnStats objects for each collision reaction contributing to k1 [NOT IMPLEMENTED]
        - ComplexRxnStats objects for each sub-reaction contributing to k2 [NOT IMPLEMENTED]
  """
  def __init__(self, *args, **kargs):
    ## Extract function arguments
    reactants = kargs['reactants']
    products = kargs['products']
    multijob = kargs.get('multistrand_job', None)
    multijob_tag = kargs.get('tag', 'success')
  
    ## Store reactants and products
    self.reactants = reactants[:]
    self.products = products[:]
      
    ## Make Multistrand job object for k1 and k2
    if multijob == None:
      self.multijob = FirstStepModeJob(self.reactants, [self.products], [multijob_tag])
    else:
      self.multijob = multijob
    self.multijob_tag = multijob_tag
    
    ## Initialize RestingSetStats list
    self.rs_stats = {rs: None for rs in self.reactants}
    
    ## Initialize reaction lists
    self.k1_rxn_stats = []
    self.k2_rxn_stats = []
    
  
  def _get_reduced_k1_stats(self, relative_error, max_sims):
    """ Returns a tuple of the net k1 value and standard error,
    calculated by taking into account temporary
    reactant depletion that reduces the effective rate of this reaction. 
    Additional simulations are run until the fractional error is
    below the given relative_error threshold. """
    fractions = [1 - rs.get_temporary_depletion(relative_error, max_sims) for rs in self.rs_stats.values() if rs != None]
    rate_fraction = reduce(lambda x,y: x*y, fractions, 1.0)
    raw_k1 = self.get_k1(relative_error / rate_fraction, max_sims)
    raw_k1_err = set.get_k1_err(max_sims = 0)
    k1 = (raw_k1*rate_fraction, raw_k1_err*rate_fraction)
    return k1

  def get_reduced_k1(self, relative_error = 0.50, max_sims = 5000, **kwargs):
    """ Returns the net k1 value, calculated by taking into account temporary
    reactant depletion that reduces the effective rate of this reaction. 
    Additional simulations are run until the fractional error is
    below the given relative_error threshold. """
    return self._get_reduced_k1_stats(relative_error, max_sims, **kwargs)[0]
  def get_k1(self, relative_error = 0.50, max_sims = 5000, **kwargs):
    """ Returns the raw k1 value calculated from Multistrand trajectory
    simulations. The relative_error is the fractional error allowed in
    the return value. Additional trials are simulated until this error
    threshold is achieved. """
    return self.get_raw_stat('k1',relative_error, max_sims, **kwargs)[0]
  def get_k2(self, relative_error = 0.50, max_sims = 5000, **kwargs):
    """ Returns the k2 folding rate on successful Multistrand trajectories. 
    Additional simulations are run until the fractional error is
    below the given relative_error threshold. """
    return self.get_raw_stat('k2', relative_error, max_sims, **kwargs)[0]
  def get_kcoll(self, relative_error = 0.50, max_sims = 5000, **kwargs):
    """ Returns the average kcoll value calculated over successful Multistrand
    trajectories. Additional simulations are run until the fractional error is
    below the given relative_error threshold. """
    return self.get_raw_stat('kcoll', relative_error, max_sims, **kwargs)[0]
  def get_prob(self, relative_error = 0.50, max_sims = 5000, **kwargs):
    """ Returns the fraction of Multistrand trajectories that ended
    in the product states given to this Stats object. This may not be
    a meaningful value. Additional simulations are run until the fractional
    error is below the given relative_error threshold. """
    return self.get_raw_stat('prob', relative_error, max_sims, **kwargs)[0]

  def get_reduced_k1_error(self, relative_error = 0.50, max_sims = 5000, **kwargs):
    """ Returns the standard error on the net k1 value """
    return self._get_reduced_k1_stats(relative_error, max_sims, **kwargs)[1]
  def get_k1_error(self, relative_error = 0.50, max_sims = 5000, **kwargs):
    """ Returns the standard error on the raw k1 value """
    return self.get_raw_stat('k1', relative_error, max_sims, **kwargs)[1]
  def get_k2_error(self, relative_error = 0.50, max_sims = 5000, **kwargs):
    """ Returns the standard error on k2 """
    return self.get_raw_stat('k2', relative_error, max_sims, **kwargs)[1]
  def get_kcoll_error(self, relative_error = 0.50, max_sims = 5000, **kwargs):
    """ Returns the standard error on kcoll """
    return self.get_raw_stat('kcoll', relative_error, max_sims, **kwargs)[1]
  def get_prob_error(self, relative_error = 0.50, max_sims = 5000, **kwargs):
    """ Returns the standard error on the probability """
    return self.get_raw_stat('prob', relative_error, max_sims, **kwargs)[1]

  def get_simulation_data(self):
    return self.multijob.get_simulation_data()

  def get_num_sims(self, tag = None):
    if tag is None:
      return self.get_multistrandjob().total_sims
    else:
      tag_id = self.get_multistrandjob().tag_id_dict[tag]
      return (self.get_simulation_data()['tags'] == tag_id).sum()
  def get_num_successful_sims(self):
    return self.get_num_sims(tag = self.get_multistrand_tag())
  def get_num_failed_sims(self):
    return self.get_num_sims() - self.get_num_successful_sims() - self.get_num_timeout_sims()
  def get_num_timeout_sims(self):
    return self.get_num_sims() - self.get_simulation_data()['valid'].sum()
    
  def get_raw_stat(self, stat, relative_error, max_sims, **kwargs):
    """ General function to reduce the error on the given statistic
    to below the given threshold and return the value and standard
    error of the statistic. """
    # Reduce error to threshold
    self.multijob.reduce_error_to(relative_error, max_sims, reaction = self.multijob_tag, stat = stat, **kwargs)
    # Calculate and return statistic
    val = self.multijob.get_statistic(self.multijob_tag, stat)
    error = self.multijob.get_statistic_error(self.multijob_tag, stat)
    return (val, error)
 
  def set_rs_stats(self, reactant, stats):
    self.rs_stats[reactant] = stats
  def get_rs_stats(self, reactant):
    return self.rs_stats[reactant]
    
  def add_k1_rxn(self, rxn_stats):
    self.k1_rxn_stats.append(rxn_stats)
  def add_k2_rxn(self, rxn_stats):
    self.k2_rxn_stats.append(rxn_stats)

  def get_multistrandjob(self):
    return self.multijob
  def get_multistrand_tag(self):
    return self.multijob_tag
  
  
    
class RestingSetStats(object):
  """ Calculates statistics for this resting set.
      The following statistics are calculated directly:
        - probability of each resting set conformation
        - probability of each nucleotide being bound correctly [NOT IMPLEMENTED]
        - getting MFE structure (or top N structures)
      The following information should be stored with this object
      for convenience in the form of other stats objects,
      using the add_XXX_rxn() functions:
        - inter-RS reactions (RestingSetRxnStats)
        - spurious reactions (RestingSetRxnStats)
  """
  
  
  def __init__(self, restingset, kinda_params = {}, nupack_params = {}):
    """ Initialize a RestingSetStats object from a DNAObjects.RestingSet
    object. The stats objects for reaction data should be added later
    with the add_XXX_rxn() functions. """
    ## Store resting set
    self.restingset = restingset
    self.strands = list(restingset.complexes)[0].strands
    self.strand_seqs = [s.sequence for s in self.strands]
    
    ## Set up NUPACK sampler (for conformation probabilities)
    self.sampler = NupackSampleJob(
        restingset,
        similarity_threshold = kinda_params.get('nupack_similarity_threshold', None),
        multiprocessing = kinda_params.get('nupack_multiprocessing', True),
        nupack_params = nupack_params
    )
    
    ## Set up MFE structures list
    self.mfe_structs = []

    ## Set up parameters for certain calculations
    self.c_max = None
  
    ## Set up reaction stat lists
    self.inter_rxns = []
    self.spurious_rxns = []
    
  def get_conformation_prob(self, complex_name, relative_error = 0.50, max_sims = 100000, **kwargs):
    """ Returns the probability and probability error
    of the given conformation based on the current number of samples.
    Use None as a complex name to get the probability of
    conformations that do not match up with any of the expected
    conformations. """
    self.sampler.reduce_error_to(relative_error, max_sims, complex_name, **kwargs)
    prob = self.sampler.get_complex_prob(complex_name)
    return prob
  def get_conformation_prob_error(self, complex_name, relative_error = 0.50, max_sims = 100000, **kwargs):
    self.sampler.reduce_error_to(relative_error, max_sims, complex_name, **kwargs)
    error = self.sampler.get_complex_prob_error(complex_name)
    return error
  def get_conformation_probs(self, relative_error = 0.50, max_sims = 100000, **kwargs):
    """ Returns the probability and probability error for all
    conformations in the resting set as a dictionary. """
    names = [c.name for c in self.restingset.complexes] + [None]
    for n in names:
      self.sampler.reduce_error_to(relative_error, max_sims, n, **kwargs)
    return {name: self.get_conformation_prob(name, max_sims = 0) for name in names}
  def get_similarity_threshold(self):
    return self.sampler.similarity_threshold
  def set_similarity_threshold(self, threshold):
    self.sampler.set_similarity_threshold(threshold)

  def get_conformation_prob_data(self, complex_name):
    return self.sampler.get_complex_prob_data(complex_name)

  def get_num_simulations(self):
    return self.get_nupackjob().total_sims
    
  def get_top_MFE_structs(self, num):
    """ Attempts to obtain the top <num> MFE structures by calling
    NUPACK's subopt executable with increasingly higher energy gaps
    until enough structures are returned. """
    return self.get_nupackjob().get_top_MFE_structs(num)
    
  def get_temporary_depletion_due_to(self, rxn, relative_error = 0.5, max_sims=500):
    binding_polynomial = 1. / (1 - self.get_temporary_depletion(relative_error, max_sims))

    other_reactant = rxn.reactants[0] if rxn.reactants[0]!=self.restingset else rxn.reactants[1]
    c_max = rxn.get_rs_stats(other_reactant).c_max
    if c_max is None or c_max == 0:  return 0
    return c_max * rxn.get_k1(max_sims = 0) / rxn.get_k2(max_sims = 0) / binding_polynomial
  def get_temporary_depletion(self, relative_error = 0.5, max_sims = 500):
    binding_polynomial = 1.0

    rxns = list(filter(lambda r: len(r.reactants)>1 and r.reactants == r.products, self.inter_rxns))
    for rxn in rxns:
      other_reactant = rxn.reactants[0] if rxn.reactants[0]!=self.restingset else rxn.reactants[1]
      c_max = rxn.get_rs_stats(other_reactant).c_max
      if c_max is None or c_max == 0:  continue

      rxn.get_k1(relative_error, max_sims = max_sims)
      rxn.get_k2(relative_error, max_sims = max_sims)
      binding_polynomial += c_max * rxn.get_k1(max_sims = 0) / rxn.get_k2(max_sims = 0)
     
    return 1 - 1/binding_polynomial

  def get_permanent_depletion_due_to(self, rxn, relative_error, max_sims):
    conc = 1.0
    for reactant in rxn.reactants:
      c_max = rxn.get_rs_stats(reactant).c_max
      if c_max == None:
        return 0.0
      conc *= c_max
    dep = conc * rxn.get_k1(relative_error, max_sims = max_sims) / self.c_max
    return dep if not math.isnan(dep) else 0.0
  def get_permanent_depletion(self, relative_error = 0.5, max_sims = 500):
    """ Returns the rate in /s at which the given resting set is depleted due to spurious reactions """
    if self.c_max == None or self.c_max == 0:
      return 0.0

    depletions = [self.get_permanent_depletion_due_to(rxn, relative_error, max_sims) for rxn in self.spurious_rxns]
    return sum(depletions)
  
  def add_inter_rxn(self, rxn):
    self.inter_rxns.append(rxn)
  def add_spurious_rxn(self, rxn):
    self.spurious_rxns.append(rxn)

  def get_nupackjob(self):
    return self.sampler
  
### NOTE: ComplexRxnStats and ComplexStats are not currently fully implemented
#
# class ComplexRxnStats(object):
#   """ Calculates statistics for this complex-level reaction.
#       The following statistics are calculated directly:
#         - rate constant (via Multistrand)
#         - rate constants of base-by-base subreactions [NOT IMPLEMENTED]
#       The following information should be stored with this object for
#       convenience in displaying and calculating certain information:
#         - ComplexStats objects for each reactant
#   """
#   def __init__(self, reactants, products, loose_products = True, loose_cutoff = None):
#     ## Store reactants and products
#     self.reactants = reactants[:]
#     self.products = products[:]
#     
#     ## If necessary, create "loose" product macrostates allowing a set amount
#     ## of deviation from the exact product complexes.
#     if loose_products:
#       if loose_cutoff == None:
#         loose_cutoff = options.kinda_params['loose_complex_similarity']
#       macrostates = [utils.similar_complex_macrostate(c, loose_cutoff) for c in products]
#     else:
#       macrostates = [utils.exact_complex_macrostate(c) for c in products]
#     self.stop_conditions = Macrostate(name = 'success',
#                                                  type = 'conjunction',
#                                                  macrostates = macrostates)
#                   
#     ## Set up Multistrand simulation job for the overall reaction
#     self.rxnjob = FirstPassageTimeModeJob(self.reactants, self.stop_conditions)
#     
#     ## Set up Multistrand simulation jobs for base-by-base reactions
#     ## [NOT IMPLEMENTED]
#     
#     ## Initialize ComplexStats list for reactants
#     self.complex_stats = dict([(c.id, None) for c in reactants])
#     
#   def get_rate(self, relative_error = 0.50, max_sims = 500):
#     self.rxnjob.reduce_error_to(relative_error, max_sims, 'success', 'rate')
#     rate = self.rxnjob.get_statistic('success', 'rate')
#     error = self.rxnjob.get_statistic_error('success', 'rate')
#     return (rate, error)
#     
#   def set_complex_stats(self, complex_id, stats):
#     self.complex_stats[complex_id] = stats
#   def get_complex_stats(self, complex_id):
#     return self.complex_stats[complex_id]
#   
# class ComplexStats(object):
#   pass
# 

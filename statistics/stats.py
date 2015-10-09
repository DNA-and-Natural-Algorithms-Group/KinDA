#### TODO: RestingSetRxnStats should be updated to take a DNAObjects.Reaction object,
####       rather than reactants and products separately.

## IMPORTS

from simulation.multistrandjob import FirstPassageTimeModeJob, \
    TransitionModeJob, FirstStepModeJob
from simulation.nupackjob import NupackSampleJob
from enumeration.enumeratejob import EnumerateJob
import stats_utils
import options
  
from dnaobjects import utils, Macrostate

#import PyNupack as nupack
from imports import PyNupack as nupack
    

## GLOBALS

## CLASSES
class SystemStats(object):
  """ Stores and manages Stats objects for each system component for easier retrieval.
      Data can also be stored in a file and retrieved for later analysis.
      A SystemStats object is instantiated with an EnumerateJob object, from which
      detailed and condensed reactions as well as resting sets and complexes are taken. """

  reactions = None
  spurious_reactions = None

  rs_reactions = None
  spurious_rs_reactions = None

  restingsets = None
  complexes = None

  enum_job = None

  rxn_to_stats = None
  rs_to_stats = None

  def __init__(self, complexes, c_max = None):
    ## 
    self.enum_job = EnumerateJob(complexes = complexes)

    self.complexes = complexes
    self.restingsets = self.enum_job.get_restingsets()

    self.reactions = self.enum_job.get_reactions()
    self.rs_reactions = self.enum_job.get_restingset_reactions()

    self.rxn_to_stats = stats_utils.make_RestingSetRxnStats(self.enum_job)
    self.rs_to_stats = stats_utils.make_RestingSetStats(self.restingsets)
    
    self.spurious_reactions = []
    self.spurious_rs_reactions = list(set(self.rxn_to_stats.keys()) - set(self.rs_reactions))

    for rxn in self.rs_reactions:
      rxn_stats = self.rxn_to_stats[rxn]
      for reactant in rxn.reactants:
        rs_stats = self.rs_to_stats[reactant]
        rxn_stats.set_rs_stats(reactant, rs_stats)
        rs_stats.add_inter_rxn(rxn_stats)

    for rxn in self.spurious_rs_reactions:
      rxn_stats = self.rxn_to_stats[rxn]
      for reactant in rxn.reactants:
        rs_stats = self.rs_to_stats[reactant]
        rxn_stats.set_rs_stats(reactant, rs_stats)
        rs_stats.add_spurious_rxn(rxn_stats)

    for rs_stats in self.rs_to_stats.values():
      rs_stats.c_max = c_max



  def get_reaction(self, reactants, products):
    """ Returns a single reaction with exactly the same reactants and products as those given. """
    rxns = self.get_reactions(reactants, products)
    rxns = filter(lambda x: x.reactants_equal(reactants) and x.products_equal(products), rxns)
    if len(rxns) >= 1:
      return next(iter(rxns))
    else:
      print "Warning: Could not find reaction with the given reactants and products"
      return None

  def get_reactions(self, reactants = [], products = [], spurious = None):
    """ Returns a list of all reactions including the given reactants and the given products.
        If specified, spurious = True will return only spurious reactions (those not enumerated by Peppercorn)
        and spurious = False will return only enumerated reactions. Otherwise, no distinction will be made.
        """
    if spurious == True:
      all_rxns = self.spurious_reactions + self.spurious_rs_reactions
    elif spurious == False:
      all_rxns = self.reactions + self.rs_reactions
    else:
      all_rxns = self.spurious_reactions + self.spurious_rs_reactions + self.reactions + self.rs_reactions

    return filter(lambda x: x.has_reactants(reactants) and x.has_products(products), all_rxns)

  def get_restingsets(self, complex = None, strands = []):
    rs = self.restingsets
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

  def export_data(self, filepath):
    """ Exports all simulated data so that it can be imported in a later Python session.
    NOT IMPLEMENTED """
    pass
  def import_data(self, filepath):
    """ NOT IMPLEMENTED """
    pass



class RestingSetRxnStats(object):
  """ Calculates statistics for this resting-set reaction.
      The following statistics are calculated directly:
        - rate (via Multistrand) of collision leading to success (k1)
        - minimum k1 rate based on worst-case reactant depletion
        - rate (via Multistrand) of successful folding after collision (k2)
      The following information should be stored with this object for
      convenience in displaying and calculating certain information:
        - RestingSetStats objects for each reactant
        - ComplexRxnStats objects for each collision reaction contributing to k1
        - ComplexRxnStats objects for each sub-reaction contributing to k2
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
    
  
  def get_reduced_k1(self, allowed_error = 0.50, max_sims = 500):
    """ Returns the net k1 value, calculated by taking into account temporary
    reactant depletion that reduces the effective rate of this reaction. 
    Additional simulations are run until the fractional error is
    below the given allowed_error threshold. """
    fractions = [1 - rs.get_temp_depletion(allowed_error) for rs in self.rs_stats.values() if rs != None]
    rate_fraction = reduce(lambda x,y: x*y, fractions, 1.0)
    raw_k1 = self.get_k1(allowed_error / rate_fraction, max_sims)
    k1 = (raw_k1[0]*rate_fraction, raw_k1[1]*rate_fraction)
    return k1
  def get_k1(self, allowed_error = 0.50, max_sims = 500):
    """ Returns the raw k1 value calculated from Multistrand trajectory
    simulations. The allowed_error is the fractional error allowed in
    the return value. Additional trials are simulated until this error
    threshold is achieved. """
    return self.get_raw_stat('k1', allowed_error, max_sims)
  def get_k2(self, allowed_error = 0.50, max_sims = 500):
    """ Returns the k2 folding rate on successful Multistrand trajectories. 
    Additional simulations are run until the fractional error is
    below the given allowed_error threshold. """
    return self.get_raw_stat('k2', allowed_error, max_sims)
    
  def get_raw_stat(self, stat, allowed_error, max_sims):
    """ General function to reduce the error on the given statistic
    to below the given threshold and return the value and standard
    error of the statistic. """
    # Reduce error to threshold
    self.multijob.reduce_error_to(allowed_error, max_sims, self.multijob_tag, stat)
    # Calculate and return statistic
    val = self.multijob.get_statistic(self.multijob_tag, stat)
    error = self.multijob.get_statistic_error(self.multijob_tag, stat)
    return (val, error)
  def get_kcoll(self, allowed_error = 0.50, max_sims = 500):
    """ Returns the average kcoll value calculated over successful Multistrand
    trajectories. Additional simulations are run until the fractional error is
    below the given allowed_error threshold. """
    return self.get_raw_stat('kcoll', allowed_error, max_sims)
  def get_prob(self, allowed_error = 0.50, max_sims = 500):
    """ Returns the fraction of Multistrand trajectories that ended
    in the product states given to this Stats object. This may not be
    a meaningful value. Additional simulations are run until the fractional
    error is below the given allowed_error threshold. """
    return self.get_raw_stat('prob', allowed_error, max_sims)
  
  def set_rs_stats(self, reactant, stats):
    self.rs_stats[reactant] = stats
  def get_rs_stats(self, reactant):
    return self.rs_stats[reactant]
    
  def add_k1_rxn(self, rxn_stats):
    self.k1_rxn_stats.append(rxn_stats)
  def add_k2_rxn(self, rxn_stats):
    self.k2_rxn_stats.append(rxn_stats)
  
  
    
class RestingSetStats(object):
  """ Calculates statistics for this resting set.
      The following statistics are calculated directly:
        - probability of each resting set conformation
        - probability of each nucleotide being bound correctly [NOT IMPLEMENTED]
        - getting MFE structure (or top N structures)
      The following information should be stored with this object
      for convenience in the form of other stats objects,
      using the add_XXX_rxn() functions:
        - intra-RS reactions (ComplexRxnStats)
        - inter-RS reactions (RestingSetRxnStats)
        - spurious reactions (RestingSetRxnStats)
  """
  
  
  def __init__(self, restingset):
    """ Initialize a RestingSetStats object from a DNAObjects.RestingSet
    object. The stats objects for reaction data should be added later
    with the add_XXX_rxn() functions. """
    ## Store resting set
    self.restingset = restingset
    self.strands = list(restingset.complexes)[0].strands
    self.strand_seqs = [s.constraints for s in self.strands]
    
    ## Set up NUPACK sampler (for conformation probabilities)
    self.sampler = NupackSampleJob(restingset)
    
    ## Set up MFE structures list
    self.mfe_structs = []

    ## Set up parameters for certain calculations
    self.c_max = None
  
    ## Set up reaction stat lists
    self.intra_rxns = []
    self.inter_rxns = []
    self.spurious_rxns = []
    
  def get_conformation_prob(self, complex_name, allowed_error = 0.50, max_sims = 5000):
    """ Returns the probability and probability error
    of the given conformation based on the current number of samples.
    Use '_spurious' as a complex name to get the probability of
    conformations that do not match up with any of the expected
    conformations. """
    self.sampler.reduce_error_to(allowed_error, max_sims, complex_name)
    prob = self.sampler.get_complex_prob(complex_name)
    error = self.sampler.get_complex_prob_error(complex_name)
    return (prob, error)
  def get_conformation_probs(self, allowed_error = 0.50, max_sims = 5000):
    """ Returns the probability and probability error for all
    conformations in the resting set as a dictionary. """
    names = [c.name for c in self.restingset.complexes] + ["_spurious"]
    for n in names:
      self.sampler.reduce_error_to(allowed_error, max_sims, n)
    return {name: self.get_conformation_prob(name, max_sims = 0) for name in names}
    
  def get_top_MFE_structs(self, num):
    """ Attempts to obtain the top <num> MFE structures by calling
    NUPACK's subopt executable with increasingly higher energy gaps
    until enough structures are returned. """
    T = options.nupack_params['temp']
    material = options.nupack_params['material']
    dangles = options.nupack_params['dangles']
    energy_gap = 0.1
    struct_list = []
    while len(struct_list) < num:
      struct_list = nupack.getSubopt(self.strand_seqs, T, material, dangles, energy_gap)
      energy_gap += 0.5
    return struct_list
    
  def get_temp_depletion_due_to(self, rxn, allowed_error, max_sims):
    t_unbound = self.c_max / rxn.get_k1(allowed_error, max_sims = max_sims)[0]
    for reactant in rxn.reactants:
      c_max = rxn.get_rs_stats(reactant).c_max
      if c_max == None or c_max == 0:
        return 0.0
      t_unbound /= c_max
    t_bound = 1.0 / rxn.get_k2(allowed_error, max_sims = max_sims)[0]
    return t_bound / (t_bound + t_unbound)
  def get_temp_depletion(self, allowed_error = 0.5, max_sims = 500):
    if self.c_max == None or self.c_max == 0:
      return 0.0
      
    rxns = list(filter(lambda r: r.reactants == r.products, self.inter_rxns))
    depletions = [self.get_temp_depletion_due_to(rxn, allowed_error, max_sims) for rxn in rxns]
    return min(1, sum(depletions)) # The sum is a crude estimate of worst-case depletion

  def get_perm_depletion_due_to(self, rxn, allowed_error, max_sims):
    conc = 1.0
    for reactant in rxn.reactants:
      c_max = rxn.get_rs_stats(reactant).c_max
      if c_max == None:
        return 0.0
      conc *= c_max
    return conc * rxn.get_k1(allowed_error, max_sims = max_sims)[0] / self.c_max
  def get_perm_depletion(self, allowed_error = 0.5, max_sims = 500):
    """ Returns the rate in /s at which the given resting set is depleted due to spurious reactions """
    if self.c_max == None or self.c_max == 0:
      return 0.0

    depletions = [self.get_perm_depletion_due_to(rxn, allowed_error, max_sims) for rxn in self.spurious_rxns]
    return sum(depletions)
  
  def add_intra_rxn(self, rxn):
    self.intra_rxns.append(rxn)
  def add_inter_rxn(self, rxn):
    self.inter_rxns.append(rxn)
  def add_spurious_rxn(self, rxn):
    self.spurious_rxns.append(rxn)
  
class ComplexRxnStats(object):
  """ Calculates statistics for this complex-level reaction.
      The following statistics are calculated directly:
        - rate constant (via Multistrand)
        - rate constants of base-by-base subreactions [NOT IMPLEMENTED]
      The following information should be stored with this object for
      convenience in displaying and calculating certain information:
        - ComplexStats objects for each reactant
  """
  def __init__(self, reactants, products, loose_products = True, loose_cutoff = None):
    ## Store reactants and products
    self.reactants = reactants[:]
    self.products = products[:]
    
    ## If necessary, create "loose" product macrostates allowing a set amount
    ## of deviation from the exact product complexes.
    if loose_products:
      if loose_cutoff == None:
        loose_cutoff = options.general_params['loose_complex_similarity']
      macrostates = [utils.similar_complex_macrostate(c, loose_cutoff) for c in products]
    else:
      macrostates = [utils.exact_complex_macrostate(c) for c in products]
    self.stop_conditions = Macrostate(name = 'success',
                                                 type = 'conjunction',
                                                 macrostates = macrostates)
                  
    ## Set up Multistrand simulation job for the overall reaction
    self.rxnjob = FirstPassageTimeModeJob(self.reactants, self.stop_conditions)
    
    ## Set up Multistrand simulation jobs for base-by-base reactions
    ## [NOT IMPLEMENTED]
    
    ## Initialize ComplexStats list for reactants
    self.complex_stats = dict([(c.id, None) for c in reactants])
    
  def get_rate(self, allowed_error = 0.50, max_sims = 500):
    self.rxnjob.reduce_error_to(allowed_error, max_sims, 'success', 'rate')
    rate = self.rxnjob.get_statistic('success', 'rate')
    error = self.rxnjob.get_statistic_error('success', 'rate')
    return (rate, error)
    
  def set_complex_stats(self, complex_id, stats):
    self.complex_stats[complex_id] = stats
  def get_complex_stats(self, complex_id):
    return self.complex_stats[complex_id]
  
class ComplexStats(object):
  pass


# stats_utils.py
# Created by Joseph Berleant, 11/11/2014
#
# Provides utility functions for making stats objects, computing system-level scores,
# and exporting/importing data from a KinDA session.

####  Future possible TODOs:
####        make_ComplexRxnStats
####        make_ComplexStats
####        calc_intended_rxn_score

import sys
import itertools as it

from .. import objects as dna
from ..objects.utils import (
  restingset_count_by_complex_macrostate, restingset_count_by_domain_macrostate)
from ..simulation.multistrandjob import FirstPassageTimeModeJob, FirstStepModeJob
from .stats import RestingSetRxnStats, RestingSetStats


class SystemStatsImportError(Exception):
  pass


## GLOBALS
def listminuslist(minuend, subtrahend):
  """
  For each element in subtrahend, removes an equal element in minuend, with
  earlier elements removed first. The order of the remaining elements in minuend
  is preserved.
  """
  difference = minuend[:]
  for elem in subtrahend:
    if elem in difference: difference.remove(elem)
  return difference
  

######################################
# Utilities for making Stats objects #
######################################

def make_RestingSetRxnStats(restingsets, detailed_rxns, condensed_rxns,
    kinda_params = {}, multistrand_params = {}):
  """
  A convenience function, creating a dict mapping reactions to stats objects
  such that all stats objects with the same reactants share a Multistrand job
  object for improved efficiency.
  """
  # print("KinDA: Constructing internal KinDA objects...\r")
  sys.stdout.flush()

  # Initialize set of spurious reactions:
  # a list of all spurious reactions possible between any reactant pair
  spurious_rxns = set([])

  # Determine all possible sets of reactants
  # Each reaction may have 1 or 2 reactants
  all_reactants = set([tuple(sorted(rs)) for rs in
                       it.product(restingsets, restingsets)])
  if kinda_params['enable_unimolecular_reactions']:
    all_reactants |= set([(r,) for r in restingsets])

  # Make a Multistrand simulation job for each reactant group
  reactants_to_mjob = {}
  for i, reactants in enumerate(all_reactants):
    # Group all products coming from these reactants together
    enum_prods = [list(rxn.products) for rxn in condensed_rxns
                  if rxn.reactants_equal(reactants)]

    # Get spurious products from these reactants
    spurious_prods = get_spurious_products(reactants, detailed_rxns, enum_prods)
    new_spurious_rxns = [
        dna.RestingSetReaction(
            reactants = reactants,
            products  = p)
        for p in spurious_prods]
    spurious_rxns |= set(new_spurious_rxns)

    # Make product tags for each product group so data can be pulled out later
    tags = [str(rxn) for rxn in condensed_rxns if rxn.reactants_equal(reactants)]
    tags += [f'_spurious({rxn!s})' for rxn in new_spurious_rxns]
    spurious_flags = [False]*len(enum_prods) + [True]*len(spurious_prods)
    
    # Make Macrostates for Multistrand stop conditions
    stop_conditions = [
      create_stop_macrostate(
        state, tag, spurious = spurious_flag, options = kinda_params)
      for state, tag, spurious_flag
      in zip(enum_prods + spurious_prods, tags, spurious_flags)
    ]
    
    # Make Boltzmann sampling selector functions for each reactant
    start_macrostate_mode = kinda_params.get('start_macrostate_mode', 'ordered-complex')
    similarity_threshold = kinda_params['multistrand_similarity_threshold']
    boltzmann_selectors = [
        create_boltzmann_selector(
            restingset, 
            mode = start_macrostate_mode,
            similarity_threshold = similarity_threshold
        )
        for restingset in reactants
    ]

    # Make Multistrand job
    multiprocessing = kinda_params.get('multistrand_multiprocessing', True)
    if len(reactants) == 2:
      job = FirstStepModeJob(
          reactants,
          stop_conditions,
          boltzmann_selectors = boltzmann_selectors,
          multiprocessing = multiprocessing,
          multistrand_params = multistrand_params
      )
    elif len(reactants) == 1:
      job = FirstPassageTimeModeJob(
          reactants,
          stop_conditions,
          unimolecular_k1_scale = kinda_params['unimolecular_k1_scale'],
          boltzmann_selectors = boltzmann_selectors,
          multiprocessing = multiprocessing,
          multistrand_params = multistrand_params
      )
    reactants_to_mjob[reactants] = job

    # print(f"KinDA: Constructing internal KinDA objects... "
    #       f"{i/len(all_reactants):%}\r")
    # sys.stdout.flush()

  # Create RestingSetRxnStats object for each reaction
  rxn_to_stats = {}
  for rxn in condensed_rxns:
    rxn_to_stats[rxn] = RestingSetRxnStats(
        reactants = rxn.reactants,
        products = rxn.products,
        multistrand_job = reactants_to_mjob[rxn.reactants],
        tag = str(rxn))

  # Create a RestingSetRxnStats object for each spurious reaction between a set
  # of reactants The "unproductive" reaction will be included either as a valid
  # reaction (if it was enumerated) or spurious if not.  However, note that by
  # default we modify Peppercorn enumeration to include all unproductive
  # reactions.
  for rxn in spurious_rxns:
    rxn_to_stats[rxn] = stats = RestingSetRxnStats(
        reactants = rxn.reactants,
        products = rxn.products,
        multistrand_job = reactants_to_mjob[rxn.reactants],
        tag = f'_spurious({rxn!s})')

  # print("KinDA: Constructing internal KinDA objects... Done!")
  return rxn_to_stats
  
def get_spurious_products(reactants, reactions, stop_states):
  """ It is desirable to have Multistrand simulations end
  when interacting complexes have deviated so much from expected
  trajectories that any calculated reaction times actually
  include undesired interactions.
  Complexes along expected reaction trajectories are split up
  in all possible ways, and expected complexes are disregarded.
  This produces all unexpected complexes produced 'one step' away
  from the expected reaction trajectories. Note that complexes 
  produced by binding of the two reactants in unexpected ways
  are classified as unproductvie unless they dissociate into an
  unenumerated strand-level complex, in which case they are
  considered to be spurious.  """

  def hashable_strand_rotation(strands):
    index = 0
    poss_starts = list(range(len(strands)))
    strands_ext = strands + strands
    while len(poss_starts) > 1 and index < len(strands):
      to_compare = [strands_ext[i + index] for i in poss_starts]
      min_strand = min(to_compare, key=lambda strand: strand.name)
      poss_starts = [i for i in poss_starts if strands_ext[i + index] == min_strand]
      index += 1      
    start = poss_starts[0]
    return tuple(strands[start:] + strands[:start])

  def hashable_state(state):
    filtered = [x for x in state if x != ()]
    return tuple(sorted([hashable_strand_rotation(f) for f in filtered]))

  def enumerate_states(init_state, reactions, enumerated):
    """ Enumerates all states reachable from init_state by following any of the
    strand-list reactions given. The results are stored in <enumerated>.
    To prevent enumerating past the expected stop conditions,
    include the stop conditions as the initial value of <enumerated>.
    Note that <enumerated> is modified in place. """
    enumerated.add(init_state)
    for rxn in reactions:
      unreacting = listminuslist(list(init_state), rxn[0])
      if len(unreacting) == len(init_state) - len(rxn[0]): # check that all reactants were present
        new_state = hashable_state(unreacting + rxn[1])
        if new_state not in enumerated:
          enumerate_states(new_state, reactions, enumerated)
    return enumerated

  def binding_spurious_states(init_state, valid_states):
    """ Return spurious states one step away from given reactants,
    produced by a binding reaction between two reactants """
    spurious_states = set([])

    # states produced through binding
    for r1 in init_state:
      for r2 in listminuslist(list(init_state), [r1]):
        rot1 = [r1[i:] + r1[:i] for i in range(len(r1))]
        rot2 = [r2[i:] + r2[:i] for i in range(len(r2))]
        unreacting = listminuslist(list(init_state), [r1, r2])
        spurious_states |= set([hashable_state(unreacting + [a+b]) for a,b in it.product(rot1, rot2)])

    spurious_states = set([s for s in spurious_states if s not in valid_states])
    return spurious_states

  def dissociation_spurious_states(init_state, valid_states):
    """ Return spurious states one step away from given reactants,
    produced by a dissocation reaction. """
    spurious_states = set([])

    # states produced through disassociation
    for r in init_state:
      unreacting = listminuslist(list(init_state), [r])
      for i in range(len(r)):
        for j in range(i, len(r)):
          new_state = hashable_state(unreacting + [r[i:j], r[j:]+r[:i]])
          spurious_states.add(new_state)

    spurious_states = set([s for s in spurious_states if s not in valid_states])
    return spurious_states

  def one_step_spurious_states(init_state, valid_states):
    return (binding_spurious_states(init_state, valid_states)
            | dissociation_spurious_states(init_state, valid_states))

  def create_restingsets(valid_objects, spurious_strands):
    strands_to_restingsets = {}

    for strands in spurious_strands:
      name = ":".join(s.name for s in strands)
      c = dna.Complex(name = "cpx_" + name, strands = strands)
      strands_to_restingsets[strands] = dna.RestingSet(name = "rs_" + name, complexes = [c])

    return strands_to_restingsets
    
  # Convert to strandlist-level objects
  strandlist_reactants = hashable_state([tuple(r.strands) for r in reactants])
  strandlist_reactions = [[[hashable_strand_rotation(r.strands) for r in rxn.reactants],
                           [hashable_strand_rotation(p.strands) for p in rxn.products]]
                          for rxn in reactions]
  strandlist_stop_states = set([hashable_state([tuple(rs.strands) for rs in state])
                            for state
                            in stop_states])
    
  # Valid states consist of all states that we explicitly do NOT wish Multistrand to halt on,
  # plus the given expected stop states.
  # This consists of those states that can be enumerated from the initial state by following
  # the given reactions and all states that can be formed from a binding reaction between
  # two reactants in the initial state.
  valid_states = enumerate_states(strandlist_reactants, strandlist_reactions, strandlist_stop_states.copy()) \
                 | binding_spurious_states(strandlist_reactants, set([]))

  # The spurious states are determined as those one step away from
  # intermediate states only.
  # Stop states are not included because once a stop state is reached,
  # the simulation should halt.
  valid_intermediates = valid_states - strandlist_stop_states
  
  # Get all spurious states that can be formed by dissociation
  # within any valid intermediate state
  spurious_states = set([])
  for state in valid_intermediates:
    spurious_states |= dissociation_spurious_states(state, valid_states)

  # Create RestingSet objects given the strands
  valid_objects = set(obj for rxn in reactions for obj in rxn.reactants + rxn.products)
  spurious_strandlists = set(obj for state in spurious_states for obj in state)
  strands_to_restingsets = create_restingsets(valid_objects, spurious_strandlists)

  # Convert back to RestingSet objects
  spurious_restingsets = [
      [strands_to_restingsets[strands] for strands in state]
      for state in spurious_states
  ]
    
  return spurious_restingsets

  
def create_stop_macrostate(state, tag, spurious, options):
  """ For most simulations, there is a specific way to produce a macrostate
  from a given system state consisting of a list of complexes/resting sets.
  This procedure is as follows:
    1. For each resting set, create a Macrostate corresponding to the value of 'stop_macrostate_mode':
       a) 'ordered-complex': Create a ORDERED-COMPLEX Macrostate corresponding to the strands of the resting set.
          Because we assume that no two resting sets share the same list of ordered strands, this
          is often sufficient. This is also the fastest mode for Multistrand simulations.
       b) 'count-by-complex': Create a DISJUNCTION macrostate of COUNT macrostates, where each
          COUNT macrostate uses a fractional cutoff of options['multistrand_similarity_threshold'].
       c) 'count-by-domain': Create a DISJUNCTION macrostate of LOOSE macrostates corresponding to
          p-approximations of each domain-level complex in the resting set.
    2. Create a CONJUNCTION macrostate corresponding to the conjunction
       of the Macrostate for each resting set in the system state.
  One way to optimize this procedure would be to allow the creation of many
  macrostates from a list of system states, where each CONJUNCTION Macrostate
  could share the underlying ORDERED-COMPLEX/LOOSE/EXACT Macrostates. It's not clear
  if this would have a significant performance increase.
  Note that there is no way to represent a macrostate consisting of states with
  2 or more of a certain complex.
  """
  mode = options['stop_macrostate_mode']

  if mode not in ['ordered-complex', 'count-by-complex', 'count-by-domain']:
    print("WARNING: Unknown stop_condition_mode '{}'. Assuming 'ordered-complex'.".format(mode))
    mode = 'ordered-complex'

  obj_to_mstate = {}
  for obj in set(state):
    assert isinstance(obj, dna.RestingSet)
    if mode == 'ordered-complex' or spurious:
      obj_to_mstate[obj] = dna.Macrostate(type = 'ordered-complex', complex = next(iter(obj.complexes)))
    elif mode == 'count-by-complex':
      defect = 1 - options['multistrand_similarity_threshold']
      obj_to_mstate[obj] = restingset_count_by_complex_macrostate(obj, defect)
    elif mode == 'count-by-domain':
      defect = 1 - options['multistrand_similarity_threshold']
      obj_to_mstate[obj] = restingset_count_by_domain_macrostate(obj, defect)
          
  macrostates = dna.Macrostate(
    name        = tag,
    type        = 'conjunction',
    macrostates = [obj_to_mstate[o] for o in state])
  return macrostates


## Boltzmann sample-and-select functions must be picklable
## to be used with the multiprocessing library
class OrderedComplexSelector:
  def __init__(self, restingset):
    self._restingset = restingset # not actually used
  def __call__(self, struct):
    return True

class CountByComplexSelector:
  def __init__(self, restingset, threshold):
    self._restingset = restingset
    self._threshold = threshold
  def __call__(self, struct):
    kinda_struct = dna.Structure(strands = self._restingset.strands, structure = struct)
    return any(
        dna.utils.defect(complex, kinda_struct) < 1-self._threshold
        for complex in self._restingset.complexes
    )

class CountByDomainSelector:
  def __init__(self, restingset, threshold):
    self._restingset = restingset
    self._threshold = threshold
  def __call__(self, struct):
    kinda_struct = dna.Structure(strands = self._restingset.strands, structure = struct)
    return any(
        dna.utils.max_domain_defect(complex, kinda_struct) < 1-self._threshold
        for complex in self._restingset.complexes
    )


def create_boltzmann_selector(restingset, mode, similarity_threshold = None):
  """ Creates a 'selector' function that takes a secondary structure
  description for this resting set and returns True if it satisfies the given
  mode and similarity threshold. Available modes are:
    ordered-complex: any secondary structure with the same ordered strands
    count-by-complex: any secondary structure where the fractional defect over
      the entire complex is within 1-similarity_threshold of any of the conformations
      in the given resting set
    count-by-domain: any secondary structure where, for at least one of the conformations
      in the resting set, the fractional defect over every domain is
      within 1-similarity_threshold 
  Each function takes a single argument, the secondary structure as a dot-paren-plus string,
  which should be taken from the space of possible secondary structures with the resting set's
  strand ordering.
  """
  if mode == 'ordered-complex':
    return OrderedComplexSelector(restingset)
  elif mode == 'count-by-complex':
    assert similarity_threshold is not None
    return CountByComplexSelector(restingset, similarity_threshold)
  elif mode == 'count-by-domain':
    assert similarity_threshold is not None
    return CountByDomainSelector(restingset, similarity_threshold)


def make_RestingSetStats(restingsets, kinda_params = {}, nupack_params = {}):
  """ A convenience function to make RestingSetStats objects for
  a list of given RestingSets. Returns a dict mapping the RestingSets
  to their corresponding stats objects. """
  rs_to_stats = {
      rs: RestingSetStats(rs, 
        kinda_params = kinda_params, 
        nupack_params = nupack_params) for rs in restingsets}
  return rs_to_stats
  

def make_stats(complexes, restingsets, detailed_rxns, condensed_rxns,
        kinda_params = {}, multistrand_params = {}, nupack_params = {}):
  """
  Creates a RestingSetRxnStats object for each resting-set reaction
  and a RestingSetStats object for each resting set.
  """

  # Make RestingSetRxnStats objects for condensed reactions and predicted
  # spurious reactions.
  rxn_to_stats = make_RestingSetRxnStats(restingsets,
          detailed_rxns, condensed_rxns, kinda_params, multistrand_params)

  # Collect all resting sets, including spurious ones predicted by make_RestingSetRxnStats()
  # and make RestingSetStats for each.
  rs_rxns = list(rxn_to_stats.keys())
  restingsets = set(restingsets + [rs for rxn in rs_rxns for rs in rxn.reactants+rxn.products])
  rs_to_stats = make_RestingSetStats(restingsets, kinda_params, nupack_params)
  
  condensed_rxns_set = set(condensed_rxns)

  # Link RestingSetStats and RestingSetRxnStats objects to each other.
  for rxn in rs_rxns:
    rxn_stats = rxn_to_stats[rxn]
    for reactant in set(rxn.reactants):
      rs_stats = rs_to_stats[reactant]
      rxn_stats.set_rs_stats(reactant, rs_stats)
      if rxn in condensed_rxns_set:
        rs_stats.add_inter_rxn(rxn_stats)
      else:
        rs_stats.add_spurious_rxn(rxn_stats)

  return rs_to_stats, rxn_to_stats


######################################
#          Score calculation         #
######################################

def calc_spurious_rxn_score(system_stats, relative_error = 0.5, max_sims = 500):
  """ Calculates the spurious reaction score, which is the maximum rate of depletion of
  any resting set due to spurious reactions. """
  max_depletion = 0.0
  for rs in system_stats.get_restingsets():
    stats = system_stats.get_stats(rs)
    max_depletion = max(max_depletion, 
            stats.get_permanent_depletion(relative_error, max_sims = max_sims))
  return max_depletion

def calc_unproductive_rxn_score(system_stats, relative_error = 0.5, max_sims = 500):
  """ Calculates the unproductive reaction score, which is the maximum fraction of any resting set
  occuped with unproductive reactions at a given time. """
  max_depletion = 0.0
  for rs in system_stats.get_restingsets():
    stats = system_stats.get_stats(rs)
    max_depletion = max(max_depletion, stats.get_temporary_depletion(relative_error, max_sims = max_sims))
  return max_depletion

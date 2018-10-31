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
import json
import pickle
import numpy as np
import itertools as it

from .. import __version__ as KINDA_VERSION
from .. import objects as dna
from .. import options
from ..simulation.multistrandjob import FirstPassageTimeModeJob, FirstStepModeJob
from .stats import RestingSetRxnStats, RestingSetStats

class SystemStatsImportError(Exception):
  pass

## GLOBALS
def listminuslist(minuend, subtrahend):
  """ For each element in subtrahend, removes an equal element in minuend, with earlier
  elements removed first. The order of the remaining elements in minuend is preserved. """
  difference = minuend[:]
  for elem in subtrahend:
    if elem in difference: difference.remove(elem)
  return difference
  

######################################
# Utilities for making Stats objects #
######################################

def make_RestingSetRxnStats(restingsets, detailed_rxns, condensed_rxns, 
    kinda_params = {}, multistrand_params = {}):
  """ A convenience function, creating a dict mapping
  reactions to stats objects such that all stats objects
  with the same reactants share a Multistrand job object
  for improved efficiency. """

  #print "KinDA: Constructing internal KinDA objects...\r",
  sys.stdout.flush()

  # Initialize set of spurious reactions
  spurious_rxns = set([]); # a list of all spurious reactions possible between any reactant pair
  
  # Determine all possible sets of reactants
  # Each reaction may have 1 or 2 reactants
  if kinda_params['enable_unimolecular_reactions']:
    all_reactants = set(
        [
            tuple(sorted([r1, r2], key = lambda rs: rs.id))
            for r1, r2
            in it.product(restingsets, restingsets)
        ] +
        [(r,) for r in restingsets]
    )
  else:
    all_reactants = set(
        [
            tuple(sorted([r1,r2], key = lambda rs: rs.id))
            for r1, r2
            in it.product(restingsets, restingsets)
        ]
    )
  
  # Make a Multistrand simulation job for each reactant group
  reactants_to_mjob = {}
  for i, reactants in enumerate(all_reactants):
    # Group all products coming from these reactants together
    enum_prods = [list(rxn.products) for rxn in condensed_rxns if rxn.reactants_equal(reactants)]

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
    tags += ['_spurious({})'.format(str(rxn)) for rxn in new_spurious_rxns]
    spurious_flags = [False]*len(enum_prods) + [True]*len(spurious_prods)
    
    # Make Macrostates for Multistrand stop conditions
    stop_conditions = [
      create_stop_macrostate(state, tag, spurious = spurious_flag, options = kinda_params)
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

    #print "KinDA: Constructing internal KinDA objects... {}%\r".format(100*i/len(all_reactants)),
    #sys.stdout.flush()
    
  # Create RestingSetRxnStats object for each reaction
  rxn_to_stats = {}
  for rxn in condensed_rxns:
    stats = RestingSetRxnStats(
        reactants = rxn.reactants,
        products = rxn.products,
        multistrand_job = reactants_to_mjob[tuple(sorted(rxn.reactants, key = lambda rs: rs.id))],
        tag = str(rxn)
    )
    rxn_to_stats[rxn] = stats

  # Create a RestingSetRxnStats object for each spurious reaction between a set
  # of reactants The "unproductive" reaction will be included either as a valid
  # reaction (if it was enumerated) or spurious if not.  However, note that by
  # default we modify Peppercorn enumeration to include all unproductive
  # reactions.
  for rxn in spurious_rxns:
    stats = RestingSetRxnStats(
        reactants = rxn.reactants,
        products = rxn.products,
        multistrand_job = reactants_to_mjob[tuple(sorted(rxn.reactants, key = lambda rs: rs.id))],
        tag = '_spurious({0})'.format(str(rxn))
    )
    rxn_to_stats[rxn] = stats

  #print "KinDA: Constructing internal KinDA objects... Done!"

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
    poss_starts = range(len(strands))
    strands_ext = strands + strands
    while len(poss_starts) > 1 and index < len(strands):
      to_compare = [strands_ext[i + index] for i in poss_starts]
      min_strand = min(to_compare, key=lambda strand: strand.name)
      poss_starts = filter(lambda i: strands_ext[i + index] == min_strand, poss_starts)
      index += 1      
    start = poss_starts[0]
    return tuple(strands[start:] + strands[:start])
  def hashable_state(state):
    filtered = filter(lambda x: x != (), state)
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
  from kinda.objects.utils import restingset_count_by_complex_macrostate, restingset_count_by_domain_macrostate

  mode = options['stop_macrostate_mode']

  if mode not in ['ordered-complex', 'count-by-complex', 'count-by-domain']:
    print "WARNING: Unknown stop_condition_mode '{}'. Assuming 'ordered-complex'.".format(mode)
    mode = 'ordered-complex'

  obj_to_mstate = {}
  for obj in set(state):
    assert(obj._object_type == 'resting-set')
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
class OrderedComplexSelector(object):
  def __init__(self, restingset):
    self._restingset = restingset # not actually used
  def __call__(self, struct):
    return True
class CountByComplexSelector(object):
  def __init__(self, restingset, threshold):
    self._restingset = restingset
    self._threshold = threshold
  def __call__(self, struct):
    kinda_struct = dna.Structure(strands = self._restingset.strands, structure = struct)
    return any(
        dna.utils.defect(complex, kinda_struct) < 1-self._threshold
        for complex in self._restingset.complexes
    )
class CountByDomainSelector(object):
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
  """ Creates a RestingSetRxnStats object for each resting-set reaction
  and a RestingSetStats object for each resting set. """

  # Make RestingSetRxnStats objects for condensed reactions and predicted spurious reactions.
  rxn_to_stats = make_RestingSetRxnStats(restingsets, 
          detailed_rxns, condensed_rxns, kinda_params, multistrand_params)

  # Collect all resting sets, including spurious ones predicted by make_RestingSetRxnStats()
  # and make RestingSetStats for each.
  rs_rxns = rxn_to_stats.keys()
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



######################################
#       Import/Export Utilities      #
######################################

def export_data(sstats, filepath, use_pickle = False):
  """ Exports data of this KinDA object so that it can be imported in a later Python session.
  Does not export the entire KinDA object (only the XXXStats data that has been collected).
  Data is exported in JSON format.
  The following constructs are exported:
    - domains
    - strands
    - complexes
    - resting sets
    - reactions
    - resting-set reactions
    - resting-set statistics objects
    - resting-set-reaction statistics objects
  Implemented, but badly and fragily dependent on KinDA object implementation """
  ## Extract all objects to be exported
  rs_reactions = sstats._condensed_reactions
  restingsets = sstats._restingsets
  reactions = sstats._detailed_reactions
  complexes = sstats._complexes
  strands = set(sum([c.strands for c in complexes], []))
  domains = set(sum([s.base_domains() for s in strands], []))

  # Note: We destroy the domain hierarchy and only store "base" domains. 
  domain_to_id = {dom: 'dom{0}_{1}'.format(dom.id, dom.name) for dom in domains} 
  # id #s are guaranteed to be unique within a Python session
  domain_to_dict = {
      d_id: {'name': d.name, 'sequence': str(d.sequence)} for d,d_id in domain_to_id.iteritems()}

  strand_to_id = {strand: 'strand{0}_{1}'.format(strand.id, strand.name) for strand in strands}
  strand_to_dict = {}
  for s, s_id in strand_to_id.iteritems():
    strand_domains = [domain_to_id[d] for d in s.base_domains()]
    strand_to_dict[s_id] = {'name': s.name, 'domains': strand_domains}

  complex_to_id = {cpx: 'cpx{0}_{1}'.format(cpx.id, cpx.name) for cpx in complexes}
  complex_to_dict = {}
  for c, c_id in complex_to_id.iteritems():
    cpx_strands = [strand_to_id[s] for s in c.strands]
    complex_to_dict[c_id] = {
        'name': c.name, 'strands': cpx_strands, 'structure': c.structure.to_dotparen()}

  rs_to_id = {rs: 'rs{0}_{1}'.format(rs.id, rs.name) for rs in restingsets}
  rs_to_dict = {}
  for rs, rs_id in rs_to_id.iteritems():
    rs_complexes = [complex_to_id[c] for c in rs.complexes]
    rs_to_dict[rs_id] = {'name': rs.name, 'complexes': rs_complexes}

  rxn_to_id = {rxn: 'rxn{0}_{1}'.format(rxn.id, str(rxn)) for rxn in reactions}
  rxn_to_dict = {}
  for r, r_id in rxn_to_id.iteritems():
    reactants = [complex_to_id[c] for c in r.reactants]
    products = [complex_to_id[c] for c in r.products]
    rxn_to_dict[r_id] = {'name': r.name, 'reactants': reactants, 'products': products}

  rsrxn_to_id = {rsrxn: 'rsrxn{0}_{1}'.format(rsrxn.id, str(rsrxn)) for rsrxn in rs_reactions}
  rsrxn_to_dict = {}
  for r, r_id in rsrxn_to_id.iteritems():
    reactants = [rs_to_id[rs] for rs in r.reactants]
    products = [rs_to_id[rs] for rs in r.products]
    rsrxn_to_dict[r_id] = {'name': r.name, 'reactants': reactants, 'products': products}

  rsstats_to_dict = {}
  for rs in restingsets:
    stats = sstats.get_stats(rs)
    rsstats_to_dict[rs_to_id[rs]] = {
        'similarity_threshold': stats.get_similarity_threshold(), 'c_max': stats.c_max}
    for c in rs.complexes:
      rsstats_to_dict[rs_to_id[rs]][complex_to_id[c]] = {
        'prob': '{0} +/- {1}'.format(
          stats.get_conformation_prob(c.name, 1, max_sims=0), 
          stats.get_conformation_prob_error(c.name, max_sims=0)),
        'similarity_data': list(stats.get_conformation_prob_data(c.name))
      }

  rsrxnstats_to_dict = {}
  for rsrxn in rs_reactions:
    stats = sstats.get_stats(rsrxn)
    sim_data = {key: d.tolist() for key,d in stats.get_simulation_data().iteritems()}
    if len(rsrxn.reactants) == 2:
      rsrxnstats_to_dict[rsrxn_to_id[rsrxn]] = {
        'prob': '{0} +/- {1}'.format(stats.get_prob(max_sims = 0), stats.get_prob_error(max_sims=0)),
        'kcoll': '{0} +/- {1}'.format(stats.get_kcoll(max_sims = 0), stats.get_kcoll_error(max_sims=0)),
        'k1': '{0} +/- {1}'.format(stats.get_k1(max_sims = 0), stats.get_k1_error(max_sims=0)),
        'k2': '{0} +/- {1}'.format(stats.get_k2(max_sims = 0), stats.get_k2_error(max_sims=0)),
        'simulation_data': sim_data,
        'invalid_simulation_data': stats.get_invalid_simulation_data(),
        'tag': stats.multijob_tag
      }
    elif len(rsrxn.reactants) == 1:
      rsrxnstats_to_dict[rsrxn_to_id[rsrxn]] = {
        'prob': '{0} +/- {1}'.format(stats.get_prob(max_sims = 0), stats.get_prob_error(max_sims=0)),
        'k1': '{0} +/- {1}'.format(stats.get_k1(max_sims = 0), stats.get_k1_error(max_sims=0)),
        'k2': '{0} +/- {1}'.format(stats.get_k2(max_sims = 0), stats.get_k2_error(max_sims=0)),
        'simulation_data': sim_data,
        'invalid_simulation_data': stats.get_invalid_simulation_data(),
        'tag': stats.multijob_tag
      }

  # Prepare the overall dict object to be JSON-ed
  sstats_dict = {
    'domains': domain_to_dict,
    'strands': strand_to_dict,
    'complexes': complex_to_dict,
    'resting-sets': rs_to_dict,
    'reactions': rxn_to_dict,
    'resting-set-reactions': rsrxn_to_dict,
    'resting-set-stats': rsstats_to_dict,
    'resting-set-reaction-stats': rsrxnstats_to_dict,
    'initialization_params': sstats.initialization_params,
    'version': KINDA_VERSION
  }
  
  if use_pickle : 
    pickle.dump(sstats_dict, open(filepath, "wb"))
  else :
    json.dump(sstats_dict, open(filepath, 'w'))

  return

def import_data(filepath, use_pickle = False):
  """ Imports a KinDA object as exported in the format specified by export_data() 

  Imports:
    - domains, strands, complexes, reactions, resting-sets, resting-set reactions
    - resting-set stats:
        => similarity-threshold 
        => c_max (concentration maximum to calculate temporary depletion)
        => list of simulation results
            - loaded int "nupackjob"
    - resting-set reaction stats:
        => load

  
  """

  if use_pickle:
    sstats_dict = pickle.load(open(filepath, "rb"))
  else:
    sstats_dict = json.load(open(filepath))

  if 'version' not in sstats_dict:
    print "# KinDA: WARNING: Imported data file has no version number. Assuming KinDA {}.".format(KINDA_VERSION)
  elif sstats_dict['version'] != KINDA_VERSION:
    print "# KinDA: WARNING: Attempting conversion from KinDA {}.".format(sstats_dict['version'], KINDA_VERSION)
    sstats_dict = _import_data_convert_version(sstats_dict, sstats_dict['version'])

  domains = {}
  for domain_id, data in sstats_dict['domains'].iteritems():
    domains[domain_id] = dna.Domain(name = data['name'], sequence = data['sequence'])

  strands = {}
  for strand_id, data in sstats_dict['strands'].iteritems():
    strand_domains = [domains[d_id] for d_id in data['domains']]
    strands[strand_id] = dna.Strand(name = data['name'], domains = strand_domains)

  complexes = {}
  for complex_id, data in sstats_dict['complexes'].iteritems():
    cpx_strands = [strands[s_id] for s_id in data['strands']]
    complexes[complex_id] = dna.Complex(name = data['name'], 
        strands = cpx_strands, structure = data['structure'])

  restingsets = {}
  for rs_id, data in sstats_dict['resting-sets'].iteritems():
    rs_complexes = [complexes[c_id] for c_id in data['complexes']]
    restingsets[rs_id] = dna.RestingSet(name = data['name'], complexes = rs_complexes)

  reactions = {}
  for rxn_id, data in sstats_dict['reactions'].iteritems():
    reactants = [complexes[c_id] for c_id in data['reactants']]
    products = [complexes[c_id] for c_id in data['products']]
    reactions[rxn_id] = dna.Reaction(name = data['name'], 
        reactants = reactants, products = products)

  rs_reactions = {}
  for rsrxn_id, data in sstats_dict['resting-set-reactions'].iteritems():
    reactants = [restingsets[rs_id] for rs_id in data['reactants']]
    products = [restingsets[rs_id] for rs_id in data['products']]
    rs_reactions[rsrxn_id] = dna.RestingSetReaction(name = data['name'], 
        reactants = reactants, products = products)
    
  kparams = sstats_dict['initialization_params']['kinda_params']
  mparams = sstats_dict['initialization_params']['multistrand_params']
  nparams = sstats_dict['initialization_params']['nupack_params']
  pparams = sstats_dict['initialization_params']['peppercorn_params']

  from .. import kinda
  sstats = kinda.System(complexes = complexes.values(), 
          restingsets = restingsets.values(), 
          detailed_reactions = reactions.values(), 
          condensed_reactions = rs_reactions.values(),
          enumeration = False,
          kinda_params = kparams, 
          multistrand_params = mparams, 
          nupack_params = nparams,
          peppercorn_params = pparams)
  
  for rs_id, data in sstats_dict['resting-set-stats'].iteritems():
    stats = sstats.get_stats(restingsets[rs_id])
    if stats == None:
      print "Warning: Could not match up stored statistics for {0}.".format(restingsets[rs_id])
      continue
    nupackjob = stats.get_nupackjob()
    num_sims = 0
    for key, val in data.iteritems():
      if key == 'similarity_threshold':
        threshold = val
      elif key == 'c_max':
        stats.c_max = val
      else:
        c = complexes[key]
        nupackjob.set_complex_prob_data(c.name, np.array(val['similarity_data']))
        #NOTE: I think it should be += here
        num_sims = len(val['similarity_data'])
    nupackjob.total_sims = num_sims
    stats.set_similarity_threshold(threshold)

  for rsrxn_id, data in sstats_dict['resting-set-reaction-stats'].iteritems():
    stats = sstats.get_stats(rs_reactions[rsrxn_id])
    if stats == None:
      print "Warning: Could not match up stored statistics for {0} with a resting-set reaction in the new KinDA object.".format(rs_reactions[rsrxn_id])
      continue
    multijob = stats.get_multistrandjob()
    stats.multijob_tag = data['tag']

    num_sims = len(data['simulation_data']['tags'])
    sim_data = {key:np.array(d) for key,d in data['simulation_data'].iteritems()}
    multijob.set_simulation_data(sim_data)
    multijob.set_invalid_simulation_data(data['invalid_simulation_data'])
    multijob.total_sims = num_sims
    
  return sstats

def _import_data_convert_version(sstats_dict, version):
  major, minor, subminor = [int(v) for v in version[1:].split('.')]
  if major == 0 and minor == 1 and subminor <= 5:
    print "KinDA: ERROR: Invalid version number {}. Conversion failed. Simulations and statistical calculations may fail.".format(version)
  elif major == 0 and minor == 1 and subminor <= 7:
    # add 'valid' entry to all simulation data
    for data in sstats_dict['resting-set-reaction-stats'].values():
      tags = data['simulation_data']['tags']
      num_sims = len(tags)
      MS_TIMEOUT, MS_ERROR = -1, -3
      data['simulation_data']['valid'] = np.array([t!=MS_TIMEOUT and t!=MS_ERROR for t in tags])
    sstats_dict['version'] = 'v0.1.10'
  elif major == 0 and minor == 1 and subminor <= 10:
    # add 'invalid_simulation_data' dict ms_results
    for data in sstats_dict['resting-set-reaction-stats'].values():
      invalid_idxs = filter(lambda i:data['simulation_data']['valid'][i]==0, range(len(data['simulation_data']['valid'])))
      data['invalid_simulation_data'] = [{'simulation_index': i} for i in invalid_idxs]
    sstats_dict['version'] = 'v0.1.11'
  elif major == 0 and minor == 1 and subminor <= 12:
    # change macrostate mode name from 'disassoc' to 'ordered-complex'
    if sstats_dict['initialization_params']['kinda_params']['start_macrostate_mode'] == 'disassoc':
      sstats_dict['initialization_params']['kinda_params']['start_macrostate_mode'] = 'ordered-complex'
    if sstats_dict['initialization_params']['kinda_params']['stop_macrostate_mode'] == 'disassoc':
      sstats_dict['initialization_params']['kinda_params']['stop_macrostate_mode'] = 'ordered-complex'
    sstats_dict['version'] = 'v0.1.13'

  return sstats_dict
    

## IMPORTS

import itertools as it

from imports import dnaobjectshome
import dnaobjects as dna
import options

from simulation.multistrandjob import FirstPassageTimeModeJob, FirstStepModeJob

from stats import RestingSetRxnStats, RestingSetStats

####  TODO: make_RestingSetStats (Is this done yet?)
####        make_ComplexRxnStats
####        make_ComplexStats
####        calc_spurious_rxn_score
####        calc_unproductive_rxn_score
####        calc_intended_rxn_score

## GLOBALS
def listminuslist(minuend, subtrahend):
  difference = minuend[:]
  for elem in subtrahend:
    if elem in difference: difference.remove(elem)
  return difference
  

######################################
# Utilities for making Stats objects #
######################################

def make_RestingSetRxnStats(enum_job):
  """ A convenience function, creating a dict mapping
  reactions to stats objects such that all stats objects
  with the same reactants share a Multistrand job object
  for improved efficiency. """
  # Pull out important items from enumeration job
  detailed_rxns = enum_job.get_reactions()
  condensed_rxns = enum_job.get_restingset_reactions()
  spurious_rxns = set([]); # a list of all spurious reactions possible between any reactant pair
  
  # Retrieve all resting sets
  restingsets = enum_job.get_restingsets()

  # Determine all possible pairs of reactants
  reactants = set([
      tuple(sorted([r1, r2], key = lambda rs: rs.id))
      for r1, r2
      in it.product(restingsets, restingsets)
  ])
  
  # Make a Multistrand simulation job for each reactant group
  reactants_to_mjob = {}
  for r in reactants:
    # Group all products coming from these reactants together
    enum_prods = [list(rxn.products) for rxn in condensed_rxns if rxn.reactants_equal(r)]

    # Get spurious products from these reactants
    spurious_prods = get_spurious_products(r, detailed_rxns, enum_prods)
    new_spurious_rxns = [
        dna.RestingSetReaction(
            reactants = r,
            products  = p)
        for p in spurious_prods]
    spurious_rxns |= set(new_spurious_rxns)
    
    # Make product tags for each product group so data can be pulled out later
    tags = [str(rxn) for rxn in condensed_rxns if rxn.reactants_equal(r)]
    tags = tags + ['_spurious({0})'.format(str(rxn)) for rxn in new_spurious_rxns]
    
    # Make Macrostates for Multistrand stop conditions
    stop_conditions = [
      create_macrostate(state, tag)
        for state, tag
        in zip(enum_prods + spurious_prods, tags)
      ]
    
    # Make Multistrand job
    job = FirstStepModeJob(r, stop_conditions)
    reactants_to_mjob[r] = job
    
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

  # Create a RestingSetRxnStats object for each spurious reaction between a set of reactants
  # The "unproductive" reaction will be included either as a valid reaction (if it was enumerated) or spurious if not
  for rxn in spurious_rxns:
    stats = RestingSetRxnStats(
        reactants = rxn.reactants,
        products = rxn.products,
        multistrand_job = reactants_to_mjob[tuple(sorted(rxn.reactants, key = lambda rs: rs.id))],
        tag = '_spurious({0})'.format(str(rxn))
    )
    rxn_to_stats[rxn] = stats

  return rxn_to_stats
  
def get_spurious_products(reactants, reactions, stop_states):
  """ It is desirable to have Multistrand simulations end
  when interacting complexes have deviated so much from expected
  trajectories that any calculated reaction times actually
  include undesired interactions.
  This is done here by calculating all ways of combining the
  (at most 2) reactants and removing ones that are expected.
  Complexes along expected reaction trajectories are split up
  in all possible ways, and expected complexes are disregarded.
  This produces all unexpected complexes produced 'one step' away
  from the expected reaction trajectories."""
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
      if len(unreacting) == len(init_state) - len(rxn[0]):
        new_state = hashable_state(unreacting + rxn[1])
        if new_state not in enumerated:
          enumerate_states(new_state, reactions, enumerated)
    return enumerated

  def binding_spurious_states(init_state, valid_states):
    # Return spurious states one step away from given reactants,
    # produced by a binding reaction between two reactants
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
    # Return spurious states one step away from given reactants,
    # produced by a dissocation reaction
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

    for obj in valid_objects:
      if obj._object_type == "complex":
        strands_to_restingsets[tuple(obj.strands)] = dna.RestingSet(complexes = [obj])
      else:
        strands_to_restingsets[tuple(c.strands)] = obj

    return strands_to_restingsets
    
    
  ## Convert to strandlist-level objects
  strandlist_reactants = hashable_state([tuple(r.strands) for r in reactants])
  strandlist_reactions = [[[hashable_strand_rotation(r.strands) for r in rxn.reactants],
                           [hashable_strand_rotation(p.strands) for p in rxn.products]]
                          for rxn in reactions]
  strandlist_stop_states = set([hashable_state([tuple(rs.strands) for rs in state])
                            for state
                            in stop_states])
    
  ## Valid states consist of all states that we explicitly do NOT wish Multistrand to halt on,
  ## plus the given expected stop states.
  ## This consists of those states those that can be enumerated from the initial state by following
  ## the given reactions and all states that can be formed from a binding reaction between
  ## two reactants in the initial state.
  valid_states = enumerate_states(strandlist_reactants, strandlist_reactions, strandlist_stop_states.copy()) \
                 | binding_spurious_states(strandlist_reactants, set([]))

  ## The spurious states are determined as those one step away from
  ## intermediate states only.
  ## Stop states are not included because once a stop state is reached,
  ## the simulation should halt.
  valid_intermediates = valid_states - strandlist_stop_states
  
  ## Get all spurious states that can be formed by dissociation
  ## within any valid intermediate state
  spurious_states = set([])
  for state in valid_intermediates:
    spurious_states |= dissociation_spurious_states(state, valid_states)
    
  # Create RestingSet objects given the strands
  valid_objects = set(sum([list(rxn.reactants + rxn.products) for rxn in reactions], []))
  spurious_strandlists = set(sum([list(state) for state in spurious_states], []))
  strands_to_restingsets = create_restingsets(valid_objects, spurious_strandlists)

  # Convert back to RestingSet objects
  spurious_restingsets = [
      [strands_to_restingsets[strands] for strands in state]
      for state in spurious_states
  ]
    
  return spurious_restingsets
  
def create_macrostate(state, tag):
  """ For most simulations, there is a specific way to produce a macrostate
  from a given system state consisting of a list of complexes/resting sets.
  This procedure is as follows:
    1. For each resting set, create a DISASSOC Macrostate corresponding to the
       strands of the resting set. This is in line with our assumption that
       no two resting sets share the same list of ordered strands.
    2. For each complex, create a LOOSE or EXACT macrostate corresonding to the
       given complex.
       If the 'loose_complex' flag is set, then a LOOSE Macrostate
       is used corresponding to a p-approximation of the domain-level complex.
       For an exact definition of p-approximation, see the paper.
       If the 'loose_complex' flag is not set, then an EXACT macrostate
       is used corresponding to the exact domain-level conformation.
    3. Create a CONJUNCTION macrostate corresponding to the conjunction
       of the Macrostate for each complex/resting set in the system state.
  One way to optimize this procedure would be to allow the creation of many
  macrostates from a list of system states, where each CONJUNCTION Macrostate
  could share the underlying DISASSOC/LOOSE/EXACT Macrostates. It's not clear
  if this would have a significant performance increase.
  Note that there is no way to represent a macrostate consisting of states with
  2 or more of a certain complex. This is a shortcoming of Multistrand
  as well as the DNAObjects package.
  """
  loose_complexes = options.flags['loose_complexes']
  loose_cutoff = options.general_params['loose_complex_similarity']
      
  obj_to_mstate = {}
  for obj in set(state):
    if obj._object_type == 'resting-set':
      obj_to_mstate[obj] = dna.Macrostate(type = 'disassoc', complex = next(iter(obj.complexes)))
    elif loose_complexes:
      obj_to_mstate[obj] = utils.similar_complex_macrostate(obj, loose_cutoff)
    else:
      obj_to_mstate[obj] = utils.exact_complex_macrostate(obj)
          
  macrostates = dna.Macrostate(
    name        = tag,
    type        = 'conjunction',
    macrostates = [obj_to_mstate[o] for o in state])
  return macrostates


def make_RestingSetStats(restingsets):
  """ A convenience function to make RestingSetStats objects for
  a list of given RestingSets. Returns a dict mapping the RestingSets
  to their corresponding stats objects. """
  rs_to_stats = {rs: RestingSetStats(rs) for rs in restingsets}
  return rs_to_stats
  
def make_stats(*args, **kargs):
  rs_rxns = kargs.get('reactions',[])
  rs = kargs.get('resting-sets', [])
  
  rxn_to_stats = make_RestingSetRxnStats(rs_rxns)
  rs_to_stats = make_RestingSetStats(rs)
  
  for rxn in rs_rxns:
    rxn_stats = rxn_to_stats[rxn]
    for reactant in rxn.reactants:
      rs_stats = rs_to_stats[reactant]
      rxn_stats.set_rs_stats(reactant, rs_stats)
      rs_stats.add_inter_rxn(rxn_stats)



######################################
#          Score calculation         #
######################################

def calc_spurious_rxn_score(system_stats, t_max, allowed_error = 0.5, max_sims = 500):
  max_depletion = 0.0
  for rs in system_stats.get_restingsets():
    stats = system_stats.get_stats(rs)
    max_depletion = max(max_depletion, min(stats.get_perm_depletion(allowed_error)*t_max, 1))
  return max_depletion

def calc_unproductive_rxn_score(system_stats, allowed_error = 0.5, max_sims = 500):
  max_depletion = 0.0
  for rs in system_stats.get_restingsets():
    stats = system_stats.get_stats(rs)
    max_depletion = max(max_depletion, stats.get_temp_depletion(allowed_error))
  return max_depletion

def calc_intended_rxn_score(system_stats):
  pass
## IMPORTS

import itertools as it

from imports import dnaobjectshome
import dnaobjects as dna
import options

from simulation.multistrandjob import FirstPassageTimeModeJob, FirstStepModeJob

from stats import RestingSetRxnStats, RestingSetStats

####  TODO: get_spurious_products is not correct. should be fixed.
####        make_RestingSetStats
####        make_ComplexRxnStats
####        make_ComplexStats

## GLOBALS
def listminuslist(minuend, subtrahend):
  difference = minuend[:]
  for elem in subtrahend:
    if elem in difference: difference.remove(elem)
  return difference
  

def make_RestingSetRxnStats(enum_job):
  """ A convenience function, creating a dict mapping
  reactions to stats objects such that all stats objects
  with the same reactants share a Multistrand job object
  for improved efficiency. """
  # Pull out important items from enumeration job
  detailed_rxns = enum_job.get_reactions()
  condensed_rxns = enum_job.get_restingset_reactions()
  
  # Determine all unique reactant groups
  reactants = set([r.reactants for r in condensed_rxns])
  
  # Make a Multistrand simulation job for each reactant group
  reactants_to_mjob = {}
  for r in reactants:
    # Group all products coming from these reactants together
    valid_prods = [list(rxn.products) for rxn in condensed_rxns if rxn.reactants == r]
    
    # Get spurious products from these reactants
    spurious_prods = get_spurious_products(r, detailed_rxns)
    
    # Make product tags for each product group so data can be pulled out later
    tags = [str(rxn) for rxn in condensed_rxns if rxn.reactants == r]
    tags = tags + ['_spurious']*len(spurious_prods)
    
    # Make Macrostates for Multistrand stop conditions
    stop_conditions = [
      create_macrostate(state, tag)
        for state, tag
        in zip(valid_prods + spurious_prods, tags)
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
        multistrand_job = reactants_to_mjob[rxn.reactants],
        tag = str(rxn)
    )
    rxn_to_stats[rxn] = stats

  # Create a RestingSetRxnStats object for each spurious reaction between a set of reactants
  for r in reactants:
    rxn = dna.RestingSetReaction(
        reactants = r,
        products = (),
        name = "_spurious"
    )
    stats = RestingSetRxnStats(
        reactants = r,
        products = (), # inelegant, but no single set of products corresponds to a spurious reaction
        multistrand_job = reactants_to_mjob[r],
        tag = "_spurious"
    )
    rxn_to_stats[rxn] = stats

  # Create a RestingSetRxnStats object for each unproductive reaction between a set of reactants
  # [NOT IMPLEMENTED] TODO

    
  return rxn_to_stats
  
def get_spurious_products(reactants, reactions):
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
  def get_valid_states(init_state, reactions, valid_states):
    valid_states.add(init_state)
    for rxn in reactions:
      unreacting = listminuslist(list(init_state), rxn[0])
      if len(unreacting) == len(init_state) - len(rxn[0]):
        new_state = hashable_state(unreacting + rxn[1])
        if new_state not in valid_states:
          get_valid_states(new_state, reactions, valid_states)
    return valid_states
  def one_step_spurious_states(init_state, valid_states):
    # Return spurious states one step away from given reactants
    one_step_states = set([])
    
    # states produced through disassociation
    for r in init_state:
      unreacting = listminuslist(list(init_state), [r])
      for i in range(len(r)):
        for j in range(i, len(r)):
          new_state = hashable_state(unreacting + [r[i:j], r[j:]+r[:i]])
          one_step_states.add(new_state)
         
    # states produced through binding
    for r1 in init_state:
      for r2 in listminuslist(list(init_state), [r1]):
        rot1 = [r1[i:] + r1[:i] for i in range(len(r1))]
        rot2 = [r2[i:] + r2[:i] for i in range(len(r2))]
        unreacting = listminuslist(list(init_state), [r1, r2])
        one_step_states |= set([hashable_state(unreacting + [a+b]) for a,b in it.product(rot1, rot2)])
    
    spurious_states = set([s for s in one_step_states if s not in valid_states])
    return spurious_states
    
    
  ## Convert to strandlist-level objects
  strandlist_reactants = hashable_state([tuple(r.strands) for r in reactants])
  strandlist_reactions = [[[hashable_strand_rotation(r.strands) for r in rxn.reactants],
                           [hashable_strand_rotation(p.strands) for p in rxn.products]]
                          for rxn in reactions]
    
  valid_states = get_valid_states(strandlist_reactants, strandlist_reactions, set([]))
  
  print valid_states
  # Get all spurious states one step away from any valid state
  spurious_states = set([])
  for s in valid_states:
    spurious_states |= one_step_spurious_states(s, valid_states)
    
  # Convert back to RestingSet objects
  spurious_restingsets = []
  for state in spurious_states:
    complexes = [dna.Complex(strands = list(s)) for s in state]
    restingsets = [dna.RestingSet(complexes = [c]) for c in complexes]
    spurious_restingsets.append(restingsets)
    
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
  at least 2 (or more) of a certain complex. This is a shortcoming of Multistrand
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
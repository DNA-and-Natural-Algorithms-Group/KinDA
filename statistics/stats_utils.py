## IMPORTS

import itertools as it

from imports import dnaobjectshome
import dnaobjects as dna

from simulation.multistrandjob import FirstPassageTimeModeJob, FirstStepModeJob

from stats import RestingSetRxnStats, RestingSetStats

####  TODO: get_spurious_products is not correct. should be fixed.
####        make_RestingSetStats
####        make_ComplexRxnStats
####        make_ComplexStats

## GLOBALS
def make_RestingSetRxnStats(reactions):
  """ A convenience function, creating a dict mapping
  reactions to stats objects such that all stats objects
  with the same reactants share a Multistrand job object
  for improved efficiency. """
  # Determine all unique reactant groups
  reactants = set([r.reactants for r in reactions])
  
  # Make a Multistrand simulation job for each reactant group
  reactants_to_mjob = {}
  for r in reactants:
    # Group all products coming from these reactants together
    valid_prods = [list(rxn.products) for rxn in reactions if rxn.reactants == r]
    
    # Get spurious products from these reactants
    spurious_prods = get_spurious_products(r, valid_prods)
    
    # Make product tags for each product group so data can be pulled out later
    tags = [str(rxn) for rxn in reactions if rxn.reactants == r]
    tags = tags + ['spurious']*len(spurious_prods)
    
    # Make Multistrand job
    job = FirstStepModeJob(r, valid_prods + spurious_prods, tags)
    reactants_to_mjob[r] = job
    
  # Create RestingSetRxnStats object for each reaction
  rxn_to_stats = {}
  for rxn in reactions:
    stats = RestingSetRxnStats(
        reactants = rxn.reactants,
        products = rxn.products,
        multistrand_job = reactants_to_mjob[rxn.reactants],
        tag = str(rxn)
    )
    rxn_to_stats[rxn] = stats
    
  return rxn_to_stats
  
  
def get_spurious_products(reactants, valid_products):
  """ It is desirable to have Multistrand simulations end
  when interacting complexes have deviated so much from expected
  trajectories that any calculated reaction times actually
  include undesired interactions. Theoretically, a set of products
  has reached this point when they must undergo one or more slow
  reactions in order to end at one of the valid product states.
  This is done here by producing all partitions of the ordered
  strand list that do not correspond to valid product states. 
  Note that this does not consider other types of slow reactions
  such as 4-way branch migration. This function also does not
  support more than 2 reactants. """
  def strand_rotations(strands):
    strands = tuple(strands)
    rotations = set()
    for i in range(len(strands)):
      rotations.add(strands[:i] + strands[i:])
    return rotations
  def strand_combinations(strands1, strands2):
    rotations1 = strand_rotations(strands1)
    rotations2 = strand_rotations(strands2)
    combos = map(lambda x: x[0] + x[1], it.product(rotations1, rotations2))
    return set(combos)
  def strand_partitions(strands):
    strands = tuple(strands)
    partitions = set()
    for i in range(len(strands)):
      for j in range(i + 1, len(strands)):
        partitions.add((strands[i:j], strands[j:] + strands[:i]))
    return partitions
  def equal_under_rotation(list1, list2):
    rotations = strand_rotations(list1)
    return any(map(lambda x: x==list2, rotations))
  def equal_under_permutation(list1, list2):
    permutations = list(it.permutations(list1))
    return any(map(lambda x: all(map(lambda y: equal_under_rotation(*y), zip(x, list2))), permutations))
  def is_spurious(strand_lists):
    # INCORRECT, talk to winfree about this
    for valid in valid_products:
      valid_strands = [rs.strands for rs in valid]
      if equal_under_permutation(valid_strands, strand_lists):
        return False
    return True
    
  if len(reactants) == 2:
    strand_lists = strand_combinations(reactants[0].strands, reactants[1].strands)
  elif len(reactants) == 1:
    strand_lists = set([tuple(reactants[0].strands)])
  else:
    raise NotImplementedError("Cannot determine spurious products of a set of more than 2 reactants.")
    
    
  partitions = set(sum(map(lambda x: list(strand_partitions(x)), strand_lists), []))
  spurious_partitions = filter(is_spurious, partitions)
  
  spurious_products = []
  for partition in spurious_partitions:
    restingsets = []
    for strandlist in partition:
      c = dna.Complex(strands = strandlist)
      rs = dna.RestingSet(complexes = [c])
      rs.name = 'spurious_' + rs.name
      restingsets.append(rs)
    spurious_products.append(restingsets)
  return spurious_products

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
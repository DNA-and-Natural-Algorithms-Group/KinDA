import sys
from typing import List

import kinda
from kinda import System
from kinda.objects import RestingSet, RestingSetReaction
from kinda.statistics.stats import RestingSetStats, RestingSetRxnStats


## Choose the accuracy target and computational budget for the rate estimation.
## Note that these options can be specified separately for each statistically
## estimated quantity.
rel_error = 0.25
max_sims = 1000


def import_system(pilpath: str) -> System:
  #### To use KinDA to analyze statistics, create a System object, which provides
  #### convenient ways to access the reactions and resting sets of a a system,
  #### along with access to the corresponding Stats objects. The System object
  #### sets a fair number of things up behind the scenes, so while it's possible
  #### to avoid using the System object (and set everything up manually), it's not
  #### recommended unless you know the code fairly well.

  #### Create the System object
  ##   Here, the System object is created automatically from a PIL file. It can
  ##   also be created with kinda.System(complexes = complexes, ...) The c_max
  ##   parameter is the default maximum concentration for any resting set, used
  ##   for calculating overall unproductive and spurious scores for the system.
  ##   c_max can be set manually for particular resting sets (see below)
  ##   [TODO: Show example]
  return kinda.from_pil(pilpath, kinda_params = {'max_concentration': 1e-7})


def export_import_data(kinda_obj: System) -> System:
  print()
  print("Exporting & importing results...")
  kinda.export_data(kinda_obj, 'analyze.db')
  kinda_obj = kinda.import_data('analyze.db')
  return kinda_obj


def get_restingsets(kinda_obj: System) -> List[RestingSet]:
  print()
  print("Resting Sets:")
  restingsets = kinda_obj.get_restingsets()
  restingsets = list(sorted(restingsets))
  for i, rs in enumerate(restingsets):
    print(f"\t{i}: {str(rs):{max(len(str(rs_)) for rs_ in restingsets)}}"
          f" = {[c.canonical_form for c in rs.complexes]}")
  return restingsets


def get_reactions(kinda_obj: System, restingsets: List[RestingSet]
                  ) -> List[RestingSetReaction]:
  ## Try:
  # Get all spurious reactions involving restingsets[0] and restingsets[1]
  rxns = kinda_obj.get_reactions(reactants=[restingsets[0], restingsets[1]],
                                 spurious=True)
  # Get all reactions
  rxns = kinda_obj.get_reactions()
  # Get only spurious reactions
  rxns = kinda_obj.get_reactions(spurious=True)

  ## Default:
  # Get all reactions that were enumerated and not unproductive
  rxns = kinda_obj.get_reactions(spurious=False, unproductive=False)
  rxns = list(sorted(rxns))

  ## Print the reactions out so you can see what's going on...
  print()
  print("Reactions between", restingsets[0], "and", restingsets[1], ":")
  for i, rxn in enumerate(rxns):
    print(f"\t{i}: {rxn}")
  return rxns


def choose_reaction(rxns: List[RestingSetReaction]) -> RestingSetReaction:
  # Select the first reaction (you can choose a different one if you like, of course)
  i_rxn = 0
  rxn = rxns[i_rxn]
  print()
  print(f"Analyzing in the following:\n\t{i_rxn}: {rxn}")
  return rxn


def sample_trajectories(rxn_stats: RestingSetRxnStats) -> None:
  print()
  print("Collecting data from Multistrand...")
  # Get an estimate for k1 with 25% error
  k1 = rxn_stats.get_k1(rel_error, max_sims=2 * max_sims, verbose=1)
  # Get an estimate for k2 with 25% error
  k2 = rxn_stats.get_k2(rel_error, max_sims=max_sims, verbose=1)
  # Estimate the probability that a random Multistrand trajectory will follow this
  # reaction (not necessarily physically significant)
  prob = rxn_stats.get_prob(rel_error, verbose=1)
  # Estimate k_coll with 25% error
  k_coll = rxn_stats.get_kcoll(rel_error, verbose=1)


def choose_restingset(restingsets: List[RestingSet]) -> RestingSet:
  i_restingset = 0
  restingset = restingsets[i_restingset]
  print()
  print(f"Analyzing in the following:\n\t{i_restingset}: {restingset}")
  ## NOTE: You can also use kinda_obj.get_restingset(strands = list_of_strands)
  ## and all resting sets involving those strands (in ANY order, and
  ## including those with additional strands) will be returned.
  return restingset


def get_conformation_probabilities(restingset: RestingSet, rs_stats: RestingSetStats):
  print()
  print(f"Conformation probabilities for resting set {restingset}:")
  # Get all conformation names in the resting set
  # (most of the time there are only 1 or 2)
  confs = [c.name for c in restingset.complexes]
  # None refers to spurious conformations
  # (that aren't similar to any expected conformations)
  confs.append(None)
  for c in confs:
    p = rs_stats.get_conformation_prob(c, .025, max_sims=max_sims)
    print(f"\t{str(c):{max(len(str(c_)) for c_ in confs)}}: {p: >7.3%}")


def get_mfe_structures(restingset: RestingSet, rs_stats: RestingSetStats):
  print()
  print(f"Top 20 MFE structures for resting set {restingset}:")
  mfe_structs = rs_stats.get_top_MFE_structs(20)
  for i,s in enumerate(mfe_structs):
    print(f"\t{i: >2d}: {s[0]} ({s[1]: 5.3f})")


def get_temporary_depletion(rs_stats: RestingSetStats):
  print()
  print("Calculating temporary depletion...")
  # Get depletion with 25% error on any relevant reaction rate
  unproductive_depletion = rs_stats.get_temporary_depletion(
    rel_error, max_sims=max_sims, verbose=1)


def get_permanent_depletion(rs_stats: RestingSetStats):
  print()
  print("Calculating permanent depletion...")
  # Get depletion with 25% error on any relevant reaction rate
  spurious_depletion = rs_stats.get_permanent_depletion(
    rel_error, max_sims=max_sims, verbose=1)

  #### To get a system-level score, use the convenience functions in stats_utils.py
  #kinda.statistics.stats_utils.calc_unproductive_rxn_score(kinda_obj)
  #kinda.statistics.stats_utils.calc_spurious_rxn_score(kinda_obj)


def main(pilpath: str):
  #### Read domains, strands, and complexes from old-style PIL file
  ## Ability to read kernel-style PIL notation to be implemented in the future.
  kinda_obj = import_system(pilpath)

  #### To analyze a reaction in detail...
  ##   1) Get the resting sets in the system (if you don't have them already)
  restingsets = get_restingsets(kinda_obj)
  ##   2) Get reactions in the system that you're interested in
  rxns = get_reactions(kinda_obj, restingsets)
  ##   3) Choose a reaction to analyze in more detail
  rxn = choose_reaction(rxns)
  ##   4) Get the associated RestingSetRxnStats object
  rxn_stats = kinda_obj.get_stats(rxn)
  ##   5) Use the Stats object to get data!
  sample_trajectories(rxn_stats)

  #### To analyze a resting set in detail...
  ##   1) Get the resting set you want
  restingset = choose_restingset(restingsets)
  ##   2) Get the associated RestingSetStats object
  rs_stats = kinda_obj.get_stats(restingset)
  ##   3) Get data!
  ## Getting the conformation probabilities
  get_conformation_probabilities(restingset, rs_stats)
  ## Getting the top 20 MFE structures
  get_mfe_structures(restingset, rs_stats)
  kinda_obj = export_import_data(kinda_obj)
  ## Getting the (fractional) reactant depletion due to unproductive reactions
  get_temporary_depletion(rs_stats)
  ## Getting the rate constant of reactant depletion due to spurious reactions (units: /s)
  get_permanent_depletion(rs_stats)
  print()
  print('Done!')


if __name__ == "__main__":
  main(sys.argv[1])

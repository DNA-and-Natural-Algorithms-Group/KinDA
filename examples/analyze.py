import sys

import kinda

#### Read domains, strands, and complexes from old-style PIL file
## Ability to read kernel-style PIL notation to be implemented in the future.
pilpath=sys.argv[1]

#### To use KinDA to analyze statistics, create a System object, which
#### provides convenient ways to access the reactions and resting sets of a
#### a system, along with access to the corresponding Stats objects.
#### The System object sets a fair number of things up behind the scenes,
#### so while it's possible to avoid using the System object (and set everything up manually),
#### it's not recommended unless you know the code fairly well.

#### Create the System object
##   Here, the System object is created automatically from a PIL file.
##   It can also be created with kinda.System(complexes = complexes, ...)
##   The c_max parameter is the default maximum concentration for any resting set, used
##   for calculating overall unproductive and spurious scores for the system.
##   c_max can be set manually for particular resting sets (see below) [TODO: Show example]
kinda_obj = kinda.from_pil(pilpath, kinda_params = {'max_concentration': 1e-7})


#### To analyze a reaction in detail...

##   1) Get the resting sets in the system (if you don't have them already)
restingsets = kinda_obj.get_restingsets() # Get all resting sets

##   2) Get reactions in the system that you're interested in
rxns = kinda_obj.get_reactions(spurious = False, unproductive=False) # Get all reactions that were enumerated and not unproductive
## Try:
# rxns = kinda_obj.get_reactions(reactants = [restingsets[0], restingsets[1]], spurious = True) # Get all spurious reactions involving restingsets[0] and restingsets[1]
# rxns = kinda_obj.get_reactions() # Get all reactions
# rxns = kinda_obj.get_reactions(spurious = True) # Get only spurious reactions

##   3) Print the reactions out so you can see what's going on...
print "Reactions between", restingsets[0], "and", restingsets[1], ":"
for i, rxn in enumerate(rxns):
  print "{0}: {1}".format(i, rxn)

##   4) Choose a reaction to analyze in more detail
rxn = rxns[0] # Select the first reaction (you can choose a different one if you like, of course)

##   5) Get the associated RestingSetRxnStats object
rxn_stats = kinda_obj.get_stats(rxn)

##   6) Use the Stats object to get data!
k1 = rxn_stats.get_k1(0.5, init_batch_size = 200) # Get an estimate for k1 with 50% error
k2 = rxn_stats.get_k2(0.25, max_sims = 500) # Get an estimate for k2 with 25% error
# prob = rxn_stats.get_prob(0.25) # Estimate the probability that a random Multistrand trajectory will follow this reaction (not necessarily physically significant)
# k_coll = rxn_stats.get_kcoll(0.25) # Estimate k_coll with 25% error


#### To analyze a resting set in detail...

##   1) Get the resting set you want
restingsets = kinda_obj.get_restingsets() # Get all resting sets
print "Resting Sets: "
for i, rs in enumerate(restingsets):
  print "{0}: {1}".format(i, rs)
restingset = restingsets[0] 
##      NOTE: You can also use kinda_obj.get_restingset(strands = list_of_strands) and all
##      resting sets involving those strands (in ANY order, and including those with additional strands)
##      will be returned.

##   2) Get the associated RestingSetStats object
rs_stats = kinda_obj.get_stats(restingset)

##   3) Get data!
## Getting the conformation probabilities
confs = [c.name for c in restingset.complexes] # Get all conformation names in the resting set (most of the time there are only 1 or 2)
confs.append(None) # None refers to spurious conformations (that aren't similar to any expected conformations)
print "Conformation probabilities for resting set {0}".format(restingset)
for c in confs:
  p = rs_stats.get_conformation_prob(c, .025)
  print "\t{0}: {1}%".format(c, p*100)
## Getting the top 10 MFE structures
mfe_structs = rs_stats.get_top_MFE_structs(10)
print "Top 10 MFE structures for resting set {0}".format(restingset)
for i,s in enumerate(mfe_structs):
  print "\t{0}: {1} ({2})".format(i, s[0], s[1])
## Getting the (fractional) reactant depletion due to unproductive reactions
# unproductive_depletion = rs_stats.get_temp_depletion(0.5) # Get depletion with 50% error on any relevant reaction rate
## Getting the rate constant of reactant depletion due to spurious reactions (units: /s)
# spurious_depletion = rs_stats.get_perm_depletion(0.5) # Get depletion with 50% error on any relevant reaction rate

#### To get a system-level score, use the convenience functions in stats_utils.py
# stats_utils.calc_unproductive_rxn_score(kinda_obj)
# stats_utils.calc_spurious_rxn_score(kinda_obj, t_max = 1.0)


# Perform data analysis for Figure TEMP (Entropy-driven catalyst, Zhang et al, Science 2007)
# -- Usage is parallel to figTEMP_simulate.py

# Usage:
#   python figTEMP_analyze.py 37 demo
#   python figTEMP_analyze.py 25 publication
#   python figTEMP_analyze.py 55 trial

exp_concs = True  # use experimental concentrations for temporary depletion, or use 100 nM default?

import sys
import datetime
import numpy
import kinda
import multistrand

# Change the MODE flag to 'demo' for quick-and-dirty (50%) results, 'trial' for more detailed (10%) results, or 'publication' for super detailed (1%) results.
MODE = 'demo'

# Read temperature
if len(sys.argv) == 1 or len(sys.argv) > 3 :
  temp = 25
else :
  if len(sys.argv) == 3 : 
    MODE = sys.argv[2]
  temp = int(sys.argv[1])

MODELIST = ['demo', 'trial', 'publication', 'randomR', 'randomT'] 
if not MODE in MODELIST:
  print 'Mode not recognized. Using demo mode.'
  MODE = 'demo'

print "Loading data for temperature T = {:d} and {} mode.".format(temp,MODE)

# Load the standard raw data file for this temperature and quality
DATA_PATH = 'figT{:d}_raw_{}.kinda'.format(temp,MODE)
ANALYSIS_PATH = 'figT{:d}_{}_analysis.csv'.format(temp,MODE)

out = open(ANALYSIS_PATH, 'w')

out.write("# Data analyzed from {} on {}\n,\n".format(DATA_PATH, datetime.datetime.now().strftime("%x at %X")))

## Import collected data.
sys_stats = kinda.import_data(DATA_PATH)

## Regardless of what default concentrations were set during the initial simulations, reset them to the desired values here.
#     Concentrations are stored as the 'c_max' field of each RestingSetStats object. By default, they are set to kinda_params['max_concentration']. 

if exp_concs:
  rs_names = ['Fuel', 'Substrate', 'Catalyst', 'Output', 'Intermediate', 'Signal', 'Waste']
  rs_c_max = [13e-9, 10e-9, 1e-9, 10e-9, 10e-9, 10e-9, 10e-9]
  for name, c_max in zip(rs_names, rs_c_max):
    sys_stats.get_stats(sys_stats.get_restingset(name = name)).c_max = c_max

## For each reaction, print k_1 and k_2 data.

## We separate the desired reactions from the unproductive reactions
rxns = sys_stats.get_reactions(spurious=False, unproductive=False)
out.write('# TARGET REACTION RATE DATA\n')
out.write('reaction,k_1,k_1 err,k_2,k_2 err,num_success,num_timeouts,num_total\n')
for r in rxns:
  rxn_stats = sys_stats.get_stats(r)

  k_1 = rxn_stats.get_k1(max_sims=0)
  k_1_err = rxn_stats.get_k1_error(max_sims=0)
  k_2 = rxn_stats.get_k2(max_sims=0)
  k_2_err = rxn_stats.get_k2_error(max_sims=0)

  num_success=rxn_stats.get_num_successful_sims()
  num_failed=rxn_stats.get_num_failed_sims()
  num_timeout=rxn_stats.get_num_timeout_sims()
  num_total=rxn_stats.get_num_sims()

  out.write('{},{},{},{},{},{},{},{}\n'.format(r, k_1, k_1_err, k_2, k_2_err,num_success,num_timeout,num_total))
out.write(',\n')

urxns = sys_stats.get_reactions(spurious=False, unproductive=True)
out.write('# UNPRODUCTIVE REACTION RATE DATA\n')
for r in [r for r in urxns if r not in rxns]:
  rxn_stats = sys_stats.get_stats(r)

  k_1 = rxn_stats.get_k1(max_sims=0)
  k_1_err = rxn_stats.get_k1_error(max_sims=0)
  k_2 = rxn_stats.get_k2(max_sims=0)
  k_2_err = rxn_stats.get_k2_error(max_sims=0)

  num_success=rxn_stats.get_num_successful_sims()
  num_failed=rxn_stats.get_num_failed_sims()
  num_timeout=rxn_stats.get_num_timeout_sims()
  num_total=rxn_stats.get_num_sims()

  out.write('{},{},{},{},{},{},{},{}\n'.format(r, k_1, k_1_err, k_2, k_2_err,num_success,num_timeout,num_total))
out.write(',\n')

## Look to see if any spurious reactions were found
srxns = sys_stats.get_reactions(spurious = True, unproductive = False)
out.write('# OBSERVED SPURIOUS REACTIONS\n')
for r in [r for r in srxns if r not in urxns]:
  rxn_stats = sys_stats.get_stats(r)

  k_1 = rxn_stats.get_k1(max_sims=0)
  k_1_err = rxn_stats.get_k1_error(max_sims=0)
  k_2 = rxn_stats.get_k2(max_sims=0)
  k_2_err = rxn_stats.get_k2_error(max_sims=0)

  num_success=rxn_stats.get_num_successful_sims()
  num_failed=rxn_stats.get_num_failed_sims()
  num_timeout=rxn_stats.get_num_timeout_sims()
  num_total=rxn_stats.get_num_sims()

  if num_success>0:
    out.write('{},{},{},{},{},{},{},{}\n'.format(r, k_1, k_1_err, k_2, k_2_err,num_success,num_timeout,num_total))
out.write(',\n')
  

## Look to see if any spurious reactions were found
srxns = sys_stats.get_reactions(spurious = True, unproductive = False)
out.write('# BOUNDED SPURIOUS REACTIONS\n')
for r in [r for r in srxns if r not in urxns]:
  rxn_stats = sys_stats.get_stats(r)

  k_1 = rxn_stats.get_k1(max_sims=0)
  k_1_err = rxn_stats.get_k1_error(max_sims=0)
  k_2 = rxn_stats.get_k2(max_sims=0)
  k_2_err = rxn_stats.get_k2_error(max_sims=0)

  num_success=rxn_stats.get_num_successful_sims()
  num_failed=rxn_stats.get_num_failed_sims()
  num_timeout=rxn_stats.get_num_timeout_sims()
  num_total=rxn_stats.get_num_sims()

  if num_success==0:
    out.write('{},{},{},{},{},{},{},{}\n'.format(r, k_1, k_1_err, k_2, k_2_err,num_success,num_timeout,num_total))
out.write(',\n')
  
rxns = sys_stats.get_reactions(spurious=False, unproductive=False)
out.write('# TARGET REACTION COMPLETION TIME HISTOGRAMS\n')
out.write('reaction,fraction1,fraction2,...,fraction30,bin0,bin1,...,bin30\n')
for r in rxns:
  rxn_stats = sys_stats.get_stats(r)
  times     = rxn_stats.get_reaction_times()
  k2        = rxn_stats.get_k2(max_sims=0)
  if len(times)==0:
    break
  freq,bins = numpy.histogram(times,bins=30,range=(0,5.0/k2))   # used to go to 1.01*times.max()
  out.write('{}'.format(r))
  for i in range(30):
    out.write(',{}'.format(freq[i]))
  for i in range(31):
    out.write(',{}'.format(bins[i]))
  out.write('\n')
out.write(',\n')


## For each resting set, print probability of expected conformations (1-p_spurious) for three NUPACK similarity scores: 0.51, 0.70, 0.90
for nupack_p in [0.51, 0.70, 0.90]:
  restingsets = sys_stats.get_restingsets()
  for rs in restingsets:
    rs_stats = sys_stats.get_stats(rs)
    rs_stats.set_similarity_threshold(nupack_p)
  restingsets = sys_stats.get_restingsets(spurious = False)

  out.write('# RESTING SET DATA WITH P = {}\n'.format(nupack_p))
  out.write('resting set,p(nonspurious),p(nonspurious) err,temporary depletion\n')
  for rs in restingsets:
    rs_stats = sys_stats.get_stats(rs)
    p_nonspurious = 1 - rs_stats.get_conformation_prob(None, max_sims=0)
    p_nonspurious_err = rs_stats.get_conformation_prob_error(None, max_sims=0)
    temp_dep = rs_stats.get_temporary_depletion(max_sims=0)
    out.write('{},{},{},{}\n'.format(rs, p_nonspurious, p_nonspurious_err, temp_dep))
  out.write(',\n')

## For each resting set, print probability of expected conformations (1-p_spurious) for all p-thresholds 0.51, 0.53, ..., 0.97, 0.99
out.write('# RESTING SET DATA WITH ALL P THRESHOLDS\n')
out.write('resting set,p(nonspurious|p=0.51),p(nonspurious|p=0.53),...,p(nonspurious|p=0.97),p(nonspurious|p=0.99)\n')

pv = [round(0.02*v + 0.51,2) for v in range(25)]

restingsets = sys_stats.get_restingsets(spurious = False)
for rs in restingsets:
  out.write('{}'.format(rs))
  for nupack_p in pv:
    rs_stats = sys_stats.get_stats(rs)
    rs_stats.set_similarity_threshold(nupack_p)
    p_nonspurious = 1 - rs_stats.get_conformation_prob(None, max_sims=0)
    out.write(',{}'.format(p_nonspurious))
  out.write('\n')
out.write(',\n')


## In a table, print out temporary depletion levels for each pairwise combination of resting sets.
out.write('# TEMPORARY DEPLETION DATA\n')
out.write((',{}'*len(restingsets)).format(*restingsets) + ',total\n')
for rs1 in restingsets:
  rs_stats = sys_stats.get_stats(rs1)
  unprod_rxns = [sys_stats.get_reactions(unproductive=True, reactants=[rs1,rs2])[0] for rs2 in restingsets]
  rxn_stats = [sys_stats.get_stats(rxn) for rxn in unprod_rxns]
  temp_deps = [rs_stats.get_temporary_depletion_due_to(s, max_sims=0) for s in rxn_stats]

  out.write(('{},'*(len(restingsets)+1)).format(rs1, *temp_deps))
  out.write('{}\n'.format(rs_stats.get_temporary_depletion(max_sims=0)))
out.write('\n')

out.close()

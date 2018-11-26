# Perform data collection for Figure 7 (Entropy-driven catalyst, Zhang et al., Science 2007)
# -- We run full system simulations, with an extra emphasis on leak, for a variety of temperatures.
# -- The user is responsible for the temperature choices, one run per temperature,
# -- and the accuracy choices, which can be demo (+- 50%), trial (+- 10%), or publication (+- 1%).
# -- Separate data files are stored for each simulation.

# Usage:
#   python figTEMP_simulate.py 37 demo
#   python figTEMP_simulate.py 25 publication
#   python figTEMP_simulate.py 55 trial

import sys,os
import kinda

# Change the MODE flag to 'demo' for quick-and-dirty (50%) results, 'trial' for more detailed (10%) results, or 'publication' for super detailed (1%) results.

if not len(sys.argv) == 3 : 
  sys.exit('ERROR: expected two arguments, temperature and simulation mode.')

temp = int(sys.argv[1])
MODE = sys.argv[2]

MODELIST = ['demo', 'trial', 'publication', 'randomR', 'randomT' ] 
if not MODE in MODELIST:
  print 'Mode not recognized. Using demo mode.'
  MODE = 'demo'

QUAL = MODE[0:3]

print "Using temperature T = {:d} and {} mode for simulations.".format(temp,MODE)

if MODE == 'demo' or MODE == 'trial' or MODE == 'publication' :
  PIL_PATH = 'Zhang2007.pil'
if MODE == 'randomR' :
  PIL_PATH = 'Zhang2007_R.pil' 
if MODE == 'randomT' :
  PIL_PATH = 'Zhang2007_T.pil' 
DATA_PATH = 'figT{:d}_raw_{}.kinda'.format(temp,MODE)

# make use of previous simulations if data is available
PREV_MODE={
  'demo':'nope','trial':'demo','publication':'trial',
  'randomR':'nope', 'randomT':'nope'
  }[MODE]
PREV_DATA_PATH = 'figT{:d}_raw_{}.kinda'.format(temp,PREV_MODE)


nupack_p = { 'dem':0.90, 'tri':0.90, 'pub':0.70, 'ran':0.70 } 
kinda_params={
  'multistrand_similarity_threshold': 0.51,      # this is the default
  'nupack_similarity_threshold': nupack_p[QUAL]  # default is 0.51
}

multistrand_params = {
  'sodium': 1.0,
  'simulation_time': 100.0, # default timeout is 1.0 seconds; some T = 15 C sims take longer than 10 sec
  'temperature': temp
}

nupack_params = {
  'sodium': 1.0,
  'T': temp
}


#### Read domains, strands, and complexes from old-style PIL file
#### Or recover that information, along with simulation results, from earlier (and perhaps lower-quality) simulations
if os.path.exists(DATA_PATH):
    print "Loading previous data for temperature T = {:d} and {} mode.".format(temp,MODE)
    sstats = kinda.import_data(DATA_PATH)
elif os.path.exists(PREV_DATA_PATH):
    print "Loading previous data for temperature T = {:d} and {} mode, boosting for {} mode.".format(temp,PREV_MODE,MODE)
    sstats = kinda.import_data(PREV_DATA_PATH)
    # The resting set analysis may have a changed analysis parameter. 
    # Since samples are saved, it is OK to update the parameter.  (Unlike for the multistrand similarity threshold.)
    restingsets = sstats.get_restingsets()
    for rs in restingsets:
      rs_stats = sstats.get_stats(rs)
      rs_stats.set_similarity_threshold(kinda_params['nupack_similarity_threshold'])
else: 
    sstats = kinda.from_pil(PIL_PATH, kinda_params=kinda_params, multistrand_params=multistrand_params, nupack_params=nupack_params)

#### Analyze reaction rates.

## Collect all reactions to simulate.
## Here we simulate all enumerated and unproductive reactions but
## do not explicitly collect data for spurious reactions, except for two known leak reactions (which are really the same)

desired_rxns = sstats.get_reactions(unproductive = False, spurious = False)

print "Desired reaction analysis order:"
for i,r in enumerate(desired_rxns):
  print i,r

rs_fuel      = sstats.get_restingset(name = 'Fuel')
rs_substrate = sstats.get_restingset(name = 'Substrate')
rs_leakint1  = sstats.get_restingset(complex_name = 'cpx_F:LB:OB', spurious=True)  # leaks Signal first
rs_leakint2  = sstats.get_restingset(complex_name = 'cpx_F:SB:LB', spurious=True)   # leaks Output first
leak_rxns    = sstats.get_reactions(reactants = [rs_fuel, rs_substrate], products = [rs_leakint1], spurious=True) \
             + sstats.get_reactions(reactants = [rs_fuel, rs_substrate], products = [rs_leakint2], spurious=True)

print "Potential leak reaction analysis order:"
for i,r in enumerate(leak_rxns):
  print i,r

unproductive_rxns = sstats.get_reactions(unproductive = True, spurious = False)

print "Unproductive reaction analysis order:"
for i,r in enumerate(unproductive_rxns):
  print i,r

print
print '==========================================='
print

######### OK, get to business now

## Set up parameters dict.
if QUAL == 'dem':
  params = {
    'relative_error':  0.5,
    'init_batch_size': 50,
    'max_batch_size':  1000,
    'max_sims':        5000,
    'sims_per_update': 1
  }
elif QUAL == 'tri':
  params = {
    'relative_error':  0.10,
    'init_batch_size': 500,
    'max_batch_size':  5000,
    'max_sims':        10000,
    'sims_per_update': 50
  }
elif QUAL == 'pub':
  params = {
    'relative_error':  0.10,
    'init_batch_size': 5000,
    'max_batch_size':  50000,
    'max_sims':        200000,
    'sims_per_update': 100
  }
elif QUAL == 'ran':
  params = {
    'relative_error':  0.25,
    'init_batch_size': 1000,
    'max_batch_size':  5000,
    'max_sims':        400000,
    'sims_per_update': 500
  }
# COMMENT on batch size: 6/4/2018 -- memory use of each spawned simulator process scales with batch size, 10000 = 1 Gig.  Better to keep batch size small, e.g. 1000 or 5000.
# With a small batch size, one can run many (e.g. 5 to 10) invocations of figTEMP_simulate.py in parallel on a 72 Gig AWS machine, without worrying about an out-of-memory crash.

# for the publication run, we have a larger max sims, but tolerate more error on k2
if QUAL == 'pub' or QUAL == 'ran':
  params_uni = dict(params)
  params_uni['relative_error']=0.50
else:
  params_uni = params

## Simulate each reaction -- desired reactions first

for i, r in enumerate(desired_rxns):
  print "\nAnalyzing desired reaction {}: {}".format(i,r)

  # Get stats object
  rxn_stats = sstats.get_stats(r)

  # Query k1 and k2 reaction rates to requested precision
  rxn_stats.get_k1(verbose=1, **params)
  rxn_stats.get_k2(verbose=1, **params_uni)
  print "k1: {} +/- {}".format(rxn_stats.get_k1(max_sims=0), rxn_stats.get_k1_error(max_sims=0))
  print "k2: {} +/- {}".format(rxn_stats.get_k2(max_sims=0), rxn_stats.get_k2_error(max_sims=0))

## Export all collected data  (we'll do it again later, but in case the system crashes...)
kinda.export_data(sstats, DATA_PATH)


## Simulate the most probable leak reactions, with the same accuracy goals.  We know they're there, so find them!

for i, r in enumerate(leak_rxns):
  print "\nAnalyzing potential leak reaction {}: {}".format(i,r)

  # Get stats object
  rxn_stats = sstats.get_stats(r)

  # Query k1 and k2 reaction rates to requested precision
  rxn_stats.get_k1(verbose=1, **params)
  rxn_stats.get_k2(verbose=1, **params_uni)
  print "k1: {} +/- {}".format(rxn_stats.get_k1(max_sims=0), rxn_stats.get_k1_error(max_sims=0))
  print "k2: {} +/- {}".format(rxn_stats.get_k2(max_sims=0), rxn_stats.get_k2_error(max_sims=0))

## Export all collected data  (we'll do it again later, but in case the system crashes...)
kinda.export_data(sstats, DATA_PATH)


## Simulate each reaction -- unproductive reactions now, with more sims, same accuracy goals

for i, r in enumerate(unproductive_rxns):
  print "\nAnalyzing unproductive reaction {}: {}".format(i,r)

  # Get stats object
  rxn_stats = sstats.get_stats(r)

  # Query k1 and k2 reaction rates to requested precision
  rxn_stats.get_k1(verbose=1, **params)
  rxn_stats.get_k2(verbose=1, **params_uni)
  print "k1: {} +/- {}".format(rxn_stats.get_k1(max_sims=0), rxn_stats.get_k1_error(max_sims=0))
  print "k2: {} +/- {}".format(rxn_stats.get_k2(max_sims=0), rxn_stats.get_k2_error(max_sims=0))

## Export all collected data  (we'll do it again later, but in case the system crashes...)
kinda.export_data(sstats, DATA_PATH)



#### Analyze resting sets.

## Collect all resting sets to analyze.
restingsets = sstats.get_restingsets()

print "Resting set analysis order:"
for i,rs in enumerate(restingsets):
  print i,rs
print

## Set up parameters dict.
if QUAL == 'dem':
  params = {
    'relative_error':  0.1,
    'init_batch_size': 50,
    'max_batch_size':  1000,
    'max_sims':        5000
  }
elif QUAL == 'tri':
  params = {
    'relative_error':  0.05,
    'init_batch_size': 100,
    'max_batch_size':  5000,
    'max_sims':        50000
  }
elif QUAL == 'pub' or QUAL == 'ran':
  params = {
    'relative_error':  0.01,
    'init_batch_size': 1000,
    'max_batch_size':  10000,
    'max_sims':        500000
  }

## Analyze each resting set
for i,rs in enumerate(restingsets):
  print "\nAnalyzing resting set {}: {}".format(i,rs)

  # Get stats object
  rs_stats = sstats.get_stats(rs)

  # Query probabilities for all conformations (excluding the spurious conformation, denoted as None)
  for c in rs.complexes:
    if not c.name == None:
      rs_stats.get_conformation_prob(c.name,verbose=1,**params)
  for c in rs.complexes:
    print c.name, rs_stats.get_conformation_prob(c.name,max_sims=0), '+/-', rs_stats.get_conformation_prob_error(c.name,max_sims=0)

## Export all collected data  (last time.  this overwrites previous exports.)
kinda.export_data(sstats, DATA_PATH)

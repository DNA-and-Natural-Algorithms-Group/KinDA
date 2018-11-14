# Perform data analysis for Figure 8 (Multiple fates system)

# Note: During Peppercorn enumeration, a warning about negative matrix entries may be reported.
#       This does not affect KinDA's analysis and can be ignored.

import kinda

# Change the MODE flag to 'demo' for quick-and-dirty results or 'publication' for more detailed results.
MODE = 'demo'

PIL_PATH = 'MultipleFates.pil'
DATA_PATH = 'fig8_unmod_raw.kinda'
# To simulate the modified system with stronger toehold GCGCGC use:
# PIL_PATH = 'MultipleFates-mod.pil'
# DATA_PATH = 'fig8_mod_raw.kinda'


## Import PIL file
print "Analyzing PIL file:", PIL_PATH
sstats = kinda.from_pil(PIL_PATH,
  kinda_params = {'max_concentration': 1e-7, 'nupack_similarity_threshold': 0.51},
  peppercorn_params = {'max_reaction_count': 2000, 'max_complex_count': 1000}
)

## Get resting sets corresponding to Gate and Interloper
rs_gate = sstats.get_restingset(complex_name = 'Gate1')
rs_interloper = sstats.get_restingset(complex_name = 'Interloper')

## Collect reactions to simulate
rxns = sstats.get_reactions(reactants = [rs_gate, rs_interloper], spurious=False, unproductive=False)

print "Reaction analysis order:"
for i,r in enumerate(rxns):
  print i,r
print

## Choose simulation parameters
if MODE == 'demo':
  params = {
    'relative_error':  0.5,
    'init_batch_size': 50,
    'min_batch_size':  50,
    'max_batch_size':  500,
    'max_sims':        500
  }
elif MODE == 'publication':
  params = {
    'relative_error':  0.05,
    'init_batch_size': 1000,
    'min_batch_size':  50,
    'max_batch_size':  1000,
    'max_sims':        2500000,
    'sims_per_worker': 10
  }

## Simulate each reaction
for i,r in enumerate(rxns):
  print "\nAnalyzing reaction {}: {}".format(i,r)

  # Get stats object
  stats = sstats.get_stats(r)

  # Query k1 and k2 reaction rates
  stats.get_k1(verbose=1, **params)
  stats.get_k2(verbose=1, **params)
  print 'k1:', stats.get_k1(max_sims=0), '+/-', stats.get_k1_error(max_sims=0)
  print 'k2:', stats.get_k2(max_sims=0), '+/-', stats.get_k2_error(max_sims=0)
  
  # Store intermediate results
  kinda.export_data(sstats, DATA_PATH)
print

## Analyze 'Gate' restingset
stats = sstats.get_stats(rs_gate)

if MODE == 'demo':
  params = {
    'relative_error': 0.5,
    'min_batch_size': 50,
    'max_batch_size': 1000,
    'max_sims':       1000
  }
elif MODE == 'publication':
  params = {
    'relative_error': 0.1,
    'min_batch_size': 50,
    'max_batch_size': 50000,
    'max_sims':       100000
  }

stats.get_conformation_probs(verbose=1, **params)

print "\nResults of analyzing 'Gate' resting set:"
cnames = [c.name for c in rs_gate.complexes] + [None]
for cname in cnames:
  p = stats.get_conformation_prob(cname, max_sims=0)
  e = stats.get_conformation_prob_error(cname, max_sims=0)
  print cname, p, '+/-', e
print

kinda.export_data(sstats, DATA_PATH)

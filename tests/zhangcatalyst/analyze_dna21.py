import os, sys

if os.path.realpath('../..') not in sys.path:
  sys.path.append(os.path.realpath('../..'))
  
from imports import dnaobjectshome

import dnaobjects as dna

from simulation import multistrandjob
from enumeration.enumeratejob import EnumerateJob
from statistics import stats_utils, stats

d, s, c = dna.io_PIL.from_PIL('./dna21_zhangcatalyst.pil')

enum_job = EnumerateJob(domains = d, strands = s, complexes = c)
rxns = enum_job.get_restingset_reactions()
restingsets = enum_job.enumerated_restingsets
  
rxn_to_stats = stats_utils.make_RestingSetRxnStats(enum_job)
rs_to_stats = stats_utils.make_RestingSetStats(restingsets)

all_rxns = sorted(rxn_to_stats.keys(), key = lambda x: [x.reactants,x.products])
rxns = sorted(rxns, key = lambda x: [x.reactants, x.products])
  

print "All Reactions: "
for i, rxn in enumerate(all_rxns):
  print "{0}: {1}".format(i, rxn)

print "Enumerated Reactions: "
for i, rxn in enumerate(rxns):
  print "{0}: {1}".format(i, rxn)
  
print "Resting Sets: "
for i, rs in enumerate(restingsets):
  print "{0}: {1}".format(i, rs)
  


#print rs_to_stats[rs[0]].get_conformation_probs()
#rxn1 = i[5][0]
#rxn2 = i[6][0]
#stats1 = rxn_to_stats[rxn1]
#stats2 = rxn_to_stats[rxn2]
#stats1.get_k1(0.05)
#stats1.get_k2(0.05)
#stats2.get_k1(0.05)
#stats2.get_k2(0.05)
#print stats1.get_k1(), stats1.get_k2()
#print stats2.get_k1(), stats2.get_k2()

sstats = stats.SystemStats(complexes = c, c_max = 1e-7)

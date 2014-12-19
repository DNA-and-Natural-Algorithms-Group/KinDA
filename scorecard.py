import __init__

import itertools as it

from imports import peppercornhome, dnaobjectshome
import enumerator.enumerator as enum
import dnaobjects as dna

from enumeration.enumeratejob import EnumerateJob
from simulation.multistrandjob import FirstStepModeJob
from statistics import stats, stats_utils

def score_PIL(filename):
  # use filename = '../KinD_scorecard/test_systems/SingleDisplace.pil'
  domains, strands, complexes = dna.io_PIL.from_PIL(filename)
  
  enum_job = EnumerateJob(domains = domains, strands = strands, complexes = complexes)
  restingset_rxns = enum_job.get_restingset_reactions()
  
  rxn_to_stats = stats_utils.make_RestingSetRxnStats(enum_job)
  
  return rxn_to_stats
  
from imports import pynupackhome, dnaobjectshome

import math

import PyNupack as nupack
from dnaobjects import utils, Complex
import options


class NupackSampleJob(object):

  restingset = None
  complexes_count = []
  total_count = 0
  
  similarity_threshold = 0.9
  
  def __init__(self, restingset, similarity_threshold = None):
    self.restingset = restingset
    self.complexes_count = [0] * (len(restingset.complexes) + 1)
    
    if similarity_threshold == None:
      self.similarity_threshold = options.general_params['loose_complex_similarity']
    else:
      self.similarity_threshold = similarity_threshold
  
  def get_complex_index(self, complex_name):
    matches = [i for i, c
                 in enumerate(self.restingset.complexes)
                 if c.name == complex_name]
    if len(matches) == 1:
      return matches[0]
    else:
      return -1
        
  def get_complex_prob(self, complex_name = "_spurious"):
    if self.total_count > 0:
      index = self.get_complex_index(complex_name)
      return float(self.complexes_count[index]) / self.total_count
    else:
      return float('nan')
        
  def get_complex_prob_error(self, complex_name = "_spurious"):
    index = self.get_complex_index(complex_name)
    prob = self.get_complex_prob(complex_name)
    if self.total_count > 0:# and prob != 0.0 and prob != 1.0:
      return math.sqrt(prob * (1 - prob) / self.total_count)
    else:
      return float('inf')
  
  def complex_is_similar(self, complex, sampled):
    return 1 - utils.max_domain_defect(sampled, complex.structure) >= self.similarity_threshold
  def sample(self, num_samples):
    T = options.nupack_params['temp']
    material = options.nupack_params['material']
    dangles = options.nupack_params['dangles']
    strands = next(iter(self.restingset.complexes)).strands
    strand_seqs = [strand.constraints
                    for strand
                    in strands]
    sampled = []
    structs = nupack.sample(num_samples, strand_seqs, T, material, dangles)
    #print structs
    for s in structs:
      sampled.append(Complex(strands = strands, structure = s))
      
    self.add_sampled_complexes(sampled)
      
  def add_sampled_complexes(self, sampled):
    for complex in sampled:
      matches = [i for i, c
                   in enumerate(self.restingset.complexes)
                   if self.complex_is_similar(c, complex)]
      for match in matches:
        self.complexes_count[match] += 1
      if len(matches) == 0:
        self.complexes_count[-1] += 1
        
      self.total_count += 1
        
    
  def reduce_error_to(self, rel_goal, abs_goal, complex_name = "_spurious"):
    error = self.get_complex_prob_error(complex_name)
    goal = max(rel_goal * self.get_complex_prob(complex_name), abs_goal)
      
    while not error <= goal:
      # Estimate additional trials based on inverse square root relationship
      # between error and number of trials
      print "Error should be reduced from %s to %s" % (error, goal)
      if error == float('inf'):
        num_trials = 250
      else:
        reduction = error / goal
        num_trials = min(250, int(self.total_count * (reduction**2 - 1) + 1))
        
      self.sample(num_trials)
      error = self.get_complex_prob_error(complex_name)
      goal = max(rel_goal * self.get_complex_prob(complex_name), abs_goal)
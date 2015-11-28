import math

#import PyNupack as nupack
from ..imports import pynupackhome, dnaobjectshome
from dnaobjects import utils, Complex

from ..imports import PyNupack as nupack
from .. import options
from datablock import Datablock


class NupackSampleJob(object):
  
  def __init__(self, restingset, similarity_threshold = None):
    self.restingset = restingset
    self.complex_tags = {tag: idx for idx, tag in enumerate([c.name for c in restingset.complexes] + ['_spurious'])}
    self.complex_counts = [0] * len(self.complex_tags)
    
    self.datablocks = {tag: Datablock() for tag in self.complex_tags if tag != '_spurious'}

    if similarity_threshold == None:  similarity_threshold = options.general_params['loose_complex_similarity']
    self.set_similarity_threshold(similarity_threshold)

    self.total_sims = 0

  def get_complex_index(self, complex_name):
    """Returns the unique index associated with this complex name (mainly for internal use). """
    assert complex_name in self.complex_tags, "Complex tag {0} not found in resting set {1}".format(complex_name, self.restingset)
    return self.complex_tags[complex_name]
        
  def get_complex_prob(self, complex_name = "_spurious"):
    index = self.get_complex_index(complex_name)
    N = self.total_sims
    Nc = self.complex_counts[index]
    return (Nc + 1.0)/(N + 2)
  def get_complex_prob_error(self, complex_name = "_spurious"):
    index = self.get_complex_index(complex_name)
    Nc = self.complex_counts[index]
    N = self.total_sims;
    return math.sqrt((Nc+1.0)*(N-Nc+1)/((N+3)*(N+2)*(N+2)));

  def get_complex_prob_data(self, complex_name = "_spurious"):
    return self.datablocks[complex_name].get_data()
  def set_complex_prob_data(self, complex_name, data):
    self.datablocks[complex_name].clear_data()
    self.datablocks[complex_name].add_data(data)

  def set_similarity_threshold(self, similarity_threshold):
    self.similarity_threshold = similarity_threshold

    ## Recalculate complex counts
    self.complex_counts = [0] * len(self.complex_tags)
    all_similarities = []
    for c in self.restingset.complexes:
      similarities = self.datablocks[c.name].get_data()
      num_similar = map(lambda v: v >= self.similarity_threshold, similarities).count(True)
      self.complex_counts[self.get_complex_index(c.name)] = num_similar
      all_similarities.append(similarities)

    ## Update complex counts
    spurious_idx = self.get_complex_index('_spurious')
    for complex_similarities in zip(*all_similarities):
      if all(map(lambda v: v < self.similarity_threshold, complex_similarities)):
        self.complex_counts[spurious_idx] += 1

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
    ## Add sampling data
    all_similarities = []
    for c in self.restingset.complexes:
      similarities = [1-utils.max_domain_defect(s, c.structure) for s in sampled]
      num_similar = map(lambda v: v >= self.similarity_threshold, similarities).count(True)
      self.complex_counts[self.get_complex_index(c.name)] += num_similar

      self.datablocks[c.name].add_data(similarities)
      all_similarities.append(similarities)

    ## Update complex counts
    spurious_idx = self.get_complex_index('_spurious')
    for complex_similarities in zip(*all_similarities):
      if all(map(lambda v: v < self.similarity_threshold, complex_similarities)):
        self.complex_counts[spurious_idx] += 1
        
    self.total_sims += len(sampled)
        
    
  def reduce_error_to(self, rel_goal, max_sims, complex_name = "_spurious"):
    num_sims = 0
    prob = self.get_complex_prob(complex_name)
    error = self.get_complex_prob_error(complex_name)
    goal = rel_goal * self.get_complex_prob(complex_name)
    while not error <= goal and num_sims < max_sims:
      # Estimate additional trials based on inverse square root relationship
      # between error and number of trials
      print "*** Conformation probability for '{0}':".format(complex_name),
      print "{0} +/- {1} (Goal: +/- {2}) [max additional sims: {3}] ***".format(prob, error, goal, max_sims - num_sims)
      if error == float('inf'):
        num_trials = 10
      else:
        reduction = error / goal
        num_trials = int(self.total_sims * (reduction**2 - 1) + 1)
        num_trials = min(500, num_trials, max_sims - num_sims, self.total_sims + 1)
        
      self.sample(num_trials)

      num_sims += num_trials
      prob = self.get_complex_prob(complex_name)
      error = self.get_complex_prob_error(complex_name)
      goal = rel_goal * self.get_complex_prob(complex_name)

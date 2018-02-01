# nupackjob.py
# Created by Joseph Berleant, 9/21/2014
#
# Implements a class for handling interactions with Nupack to obtain statistics about resting sets.


## IMPORTS

import sys
import math
import numpy as np

import multiprocessing, signal

from ..objects import utils, Complex
from .. import nupack, options
from sim_utils import print_progress_table


## GLOBALS
def sample_global(args):
  """ Global function for calling NUPACK, used for multiprocessing  """
  self = args[0]
  num_samples = args[1]

  ## Set up arguments/options for Nupack call
  strands = next(iter(self.restingset.complexes)).strands
  strand_seqs = [strand.sequence
                  for strand
                  in strands]

  ## Call Nupack
  structs = nupack.sample(num_samples, strand_seqs, **self._nupack_params)

  ## Convert each Nupack sampled structure (a dot-paren string) into a DNAObjects Complex object and process.
  sampled = [Complex(strands = strands, structure = s) for s in structs]
  return sampled



## CLASSES

class NupackSampleJob(object):
  """ The NupackSampleJob class implements basic statistics collection using Nupack.
  This class collects statistics on a particular resting set, including the probabilities
  of the constituent complexes.
  Use sample() to request a certain number of secondary structures from the Boltzmann distribution (using Nupack).
  Use get_complex_prob() to request the probability estimate for a particular complex.
  The similarity threshold may be changed with set_similarity_threshold().
  """

  def __init__(self, restingset, similarity_threshold = None, multiprocessing = True, nupack_params = {}):
    """ Constructs a NupackSampleJob object with the given resting set and similarity threshold.
    If similarity threshold is not given, the value in options.py (kinda_params['loose_complex_similarity'])
    is used. """

    # Store options
    self._multiprocessing = multiprocessing

    # Store nupack params
    self._nupack_params = dict(nupack_params)

    # Store resting set and relevant complex information
    self._restingset = restingset
    self._complex_tags = {tag: idx for idx, tag in enumerate([c.name for c in restingset.complexes] + [None])}
    self._complex_counts = [0] * len(self._complex_tags)
    
    # Data structure for storing raw sampling data
    self._data = {tag: np.array([]) for tag in self._complex_tags if tag is not None}

    self.total_sims = 0

    # Set similarity threshold, using default value in options.py if none specified
    if similarity_threshold == None:  similarity_threshold = options.kinda_params['loose_complex_similarity']
    self.set_similarity_threshold(similarity_threshold)

  @property
  def restingset(self):
    return self._restingset
  @property
  def complex_names(self):
    return sorted([c.name for c in restingset.complexes], key = lambda k: self._complex_tags[k])
  @property
  def complex_counts(self):
    return self._complex_counts[:]

  def get_complex_index(self, complex_name):
    """ Returns the unique index associated with this complex name (mainly for internal use). """
    assert complex_name in self._complex_tags, "Complex tag {0} not found in resting set {1}".format(complex_name, self.restingset)
    return self._complex_tags[complex_name]
        
  def get_complex_prob(self, complex_name = None):
    """ Returns the conformation probability associated with the given complex_name.
    complex_name should match the name field of one of the complexes in the original resting set used
    to instantiate this NupackSampleJob.
    The Bayesian estimator (n+1)/(N+2) is used to estimate this probability, in order to prevent
    improbable estimates when n=0,N.
    """
    index = self.get_complex_index(complex_name)
    N = self.total_sims
    Nc = self._complex_counts[index]
    return (Nc + 1.0)/(N + 2)
  def get_complex_prob_error(self, complex_name = None):
    """ Returns the standard error on the conformation probability associated with the given complex_name.
    complex_name should match the name field of one of the complex in the resting set used to
    construct this NupackSampleJob.
    The Bayesian estimator for the standard error is used
      SE = (n+1)*(N-n+1)/((N+3)*(N+2)*(N+2))
    to avoid artifacts when n=0,N.
    """
    index = self.get_complex_index(complex_name)
    Nc = self._complex_counts[index]
    N = self.total_sims;
    return math.sqrt((Nc+1.0)*(N-Nc+1)/((N+3)*(N+2)*(N+2)));

  def get_complex_prob_data(self, complex_name = None):
    """ Returns raw sampling data for the given complex_name.
    Data is returned as a list consisting of float values, where each float value is 1 minus the
    maximum fractional defect for any domain in that sampled secondary structure. """
    return self._data[complex_name]
  def set_complex_prob_data(self, complex_name, data):
    """ Set the raw sampling data for the given complex_name.
    Should be used only when importing an old KinDA session to restore state. """
    self._data[complex_name] = np.array(data)

  def get_num_sims(self):
    """ Returns the total number of sampled secondary structures. """
    return self.total_sims

  def set_similarity_threshold(self, similarity_threshold):
    """ Modifies the similarity threshold used to classify each sampled secondary structure as one of
    the expected complex conformations. The complex probability estimates are recomputed based
    on the new similarity threshold, so that future calls to get_complex_prob() will give
    correct statistics based on previously sampled secondary structures. """

    ## Change internal similarity threshold
    self.similarity_threshold = similarity_threshold

    ## Recalculate complex counts
    self._complex_counts = [0] * len(self._complex_tags)
    spurious_similarities = np.full(self.total_sims, True)
    for c in self.restingset.complexes:
      similarities = self._data[c.name]
      num_similar = np.sum(similarities >= self.similarity_threshold)
      self._complex_counts[self.get_complex_index(c.name)] = num_similar
      spurious_similarities *= similarities < self.similarity_threshold

    ## Update complex counts
    spurious_idx = self.get_complex_index(None)
    self._complex_counts[spurious_idx] = np.sum(spurious_similarities)

  def sample(self, **params):
    """ Calls sample_multiprocessing or sample_singleprocessing depending on the value of
    self._multiprocessing """
    if self._multiprocessing:
      self.sample_multiprocessing(**params)
    else:
      self.sample_singleprocessing(**params)

  def sample_multiprocessing(self, num_samples, samples_per_worker = None, status_func = lambda:None):
    """ Runs sample() in multiple processes. """

    # Temporarily remove SIGINT event handler from current process
    sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)

    # Create Pool instance and worker processes
    k = multiprocessing.cpu_count()
    p = multiprocessing.Pool(processes = k)

    # Restore original SIGINT event handler (if possible)
    if sigint_handler is None: sigint_handler = signal.SIG_DFL
    signal.signal(signal.SIGINT, sigint_handler)

    # Setup args for each process
    if samples_per_worker is None:
      samples_per_worker = num_samples / k
      args = [(self, samples_per_worker+1)] * (num_samples % k)
      args += [(self, samples_per_worker)] * (k - (num_samples % k))
    else:
      args = [(self, samples_per_worker)] * (num_samples / samples_per_worker)
      if num_samples%samples_per_worker > 0:  args+= [(self, num_samples%samples_per_worker)]

    it = p.imap_unordered(sample_global, args)
    p.close()

    try:
      sims_completed = 0
      for res in it:
        self.add_sampled_complexes(res)
        sims_completed += len(res)
        status_func(sims_completed)
    except KeyboardInterrupt:
      print "SIGINT: Ending NUPACK sampling prematurely..."
      p.terminate()
      p.join()
      raise KeyboardInterrupt


  def sample_singleprocessing(self, num_samples, status_func=lambda:None):
    """ Queries Nupack for num_samples secondary structures, sampled from the Boltzmann distribution
    of secondary structures for this resting set. The sampled secondary structures are automatically
    processed to compute new estimates for each conformation probability.
    The nupack_params dict that was given to this job during initialization is passed along to the
    Nupack Python interface.
    """
    results = sample_global((self, num_samples))
    self.add_sampled_complexes(results)

    status_func(len(results))

  def add_sampled_complexes(self, sampled):
    """ Processes a list of sampled Complex objects, computing the similarity to each of the
    resting set conformations and updating the estimates for each conformation probability. """

    ## For each complex in the resting set, add the number of sampled secondary structures
    ## that satisfy the similarity threshold to the running complex count.
    spurious_similarities = np.full(len(sampled), True)
    for c in self.restingset.complexes:
      similarities = np.array([1-utils.max_domain_defect(s, c.structure) for s in sampled])
      num_similar = np.sum(similarities >= self.similarity_threshold)
      self._complex_counts[self.get_complex_index(c.name)] += num_similar

      self._data[c.name] = np.concatenate((self._data[c.name], similarities))
      spurious_similarities *= similarities < self.similarity_threshold

    ## Update complex count for the spurious complex
    ## A conformation is considered spurious if it does not satisfy the similarity threshold
    ## for any of the predicted conformations in the resting set.
    spurious_idx = self.get_complex_index(None)
    self._complex_counts[spurious_idx] += np.sum(spurious_similarities)
        
    ## Update running count of number of samples
    self.total_sims += len(sampled)
        
    
  def reduce_error_to(self, rel_goal, max_sims, complex_name = None, init_batch_size = 50, min_batch_size = 50, max_batch_size = 1000):
    """ Continue querying Nupack for secondary structures sampled from the Boltzmann distribution
    until the standard error of the estimated probability for the given complex_name
    is at most rel_goal*complex_prob, or max_sims conformations have been sampled.
    If no complex_name is given, the halting condition is based on the error
    for the spurious conformation probability.
    """
    def status_func(batch_sims_done):
      prob = self.get_complex_prob(complex_name)
      error = self.get_complex_prob_error(complex_name)
      goal = rel_goal * prob
      if exp_add_sims is None:
        update_func([complex_name, prob, error, goal, "", "%d/--"%(num_sims+batch_sims_done), "--"])
      else:
        update_func([complex_name, prob, error, goal, "", "%d/%d"%(num_sims + batch_sims_done,num_sims+exp_add_sims), str(100*(num_sims+batch_sims_done)/(num_sims+exp_add_sims))+'%']) 
      
    if self._multiprocessing:
      print '[MULTIPROCESSING ON] (over %d cores)'%multiprocessing.cpu_count()
    else:
      print '[MULTIPROCESSING OFF]'


    # Get initial values
    num_sims = 0
    prob = self.get_complex_prob(complex_name)
    error = self.get_complex_prob_error(complex_name)
    goal = rel_goal * prob

    # Check if simulations are needed, return if not
    if error <= goal or num_sims >= max_sims:  return

    # Prepare progress table
    update_func = print_progress_table(
        ["complex", "prob", "error", "err goal", "", "sims done", "progress"],
        [10, 10, 10, 10, 10, 17, 10])
    update_func([complex_name, prob, error, goal, "", "--/--", "--"])

    # Run simulations
    while not error <= goal and num_sims < max_sims:
      # Estimate additional trials based on inverse square root relationship
      # between error and number of trials
      if error == float('inf'):
        num_trials = init_batch_size
        exp_add_sims = None
      else:
        reduction = error / goal
        exp_add_sims = int(self.total_sims * (reduction**2 - 1) + 1)
        num_trials = max(min(max_batch_size, exp_add_sims, max_sims - num_sims, self.total_sims + 1), min_batch_size)
        
      # Query Nupack
      self.sample(num_samples = num_trials, status_func = status_func)

      # Update estimates and goal
      num_sims += num_trials
      prob = self.get_complex_prob(complex_name)
      error = self.get_complex_prob_error(complex_name)
      goal = rel_goal * prob

    print

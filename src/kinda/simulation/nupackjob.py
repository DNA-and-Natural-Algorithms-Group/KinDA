# nupackjob.py
# Created by Joseph Berleant, 9/21/2014
#
# Implements a class for handling interactions with Nupack to obtain statistics
# about resting sets.


import math
import signal

import numpy as np
import multiprocess

import multistrand.utils.thermo as nupack

from .. import options
from ..objects import utils, Complex
from .sim_utils import print_progress_table


# NUPACK interface
def sample_global(job_spec):
  """
  Global function for calling NUPACK, used for multiprocessing.
  """
  (self, num_samples) = job_spec
  strands = next(iter(self.restingset.complexes)).strands
  strand_seqs = [strand.sequence for strand in strands]

  # Call Multistrand's Nupack wrapper
  structs = nupack.sample(strand_seqs, num_samples, **self._nupack_params)

  # Convert each Nupack sampled structure (a dot-paren string) into a
  # DNAObjects Complex object and process.
  return [Complex(strands=strands, structure=s.dp()) for s in structs]


class NupackSampleJob:
  """Calculate complex probabilities within resting sets using NUPACK.

  Args:
    restingset (dna.RestingSet()): A set of complexes with the same strand order.
    similarity_threshold: Overwrites the nupack_similarity_threshold paramter.
        One can update the similarity threshold at any time using
        set_similarity_threshold().
    multiprocessing (bool, optional): Distribute computation to all available
        cores. Defaults to True.
    nupack_params (dict): A dictionary with parameter for NUPACK.

  Use sample() to request a certain number of secondary structures from the
  Boltzmann distribution (using Nupack). Use get_complex_prob() to request the
  probability estimate for a particular complex. To update results for a new
  similarity threshold use set_similarity_threshold().
  """
  verbose = 1

  def __init__(self, restingset, similarity_threshold = None, 
               multiprocessing = True, nupack_params = {}):

    # Store options
    self.multiprocessing = multiprocessing
    self._mp_ctx = multiprocess.get_context('spawn')

    # Store nupack params
    self._nupack_params = dict(nupack_params)

    # Store resting set and relevant complex information
    self._restingset = restingset
    self._complex_tags = {tag: idx for idx, tag in enumerate(
      [c.name for c in restingset.complexes] + [None])}
    self._complex_counts = [0] * len(self._complex_tags)
    
    # Data structure for storing raw sampling data
    self._data = {tag: np.array([]) for tag in self._complex_tags if tag is not None}
    self.total_sims = 0

    # Set similarity threshold, using default value in options.py if none specified
    if similarity_threshold is None:
      similarity_threshold = options.kinda_params['nupack_similarity_threshold']
    self.set_similarity_threshold(similarity_threshold)

  @property
  def restingset(self):
    return self._restingset

  @property
  def complex_names(self):
    return sorted([c.name for c in self.restingset.complexes], key = lambda k: self._complex_tags[k])

  @property
  def complex_counts(self):
    return self._complex_counts[:]

  def get_complex_index(self, complex_name):
    """
    Returns the unique index associated with this complex name.
    """
    assert complex_name in self._complex_tags, \
      f"Complex tag {complex_name} not found in resting set {self.restingset}"
    return self._complex_tags[complex_name]
        
  def get_complex_prob(self, complex_name = None):
    """
    Returns the conformation probability associated with the given complex_name.
    complex_name should match the name field of one of the complexes in the
    original resting set used to instantiate this NupackSampleJob.
    The Bayesian estimator (n+1)/(N+2) is used to estimate this probability, in
    order to prevent improbable estimates when n=0,N.
    """
    N = self.total_sims
    Nc = self.get_complex_count(complex_name)
    return (Nc + 1.0)/(N + 2)

  def get_complex_prob_error(self, complex_name = None):
    """
    Returns the standard error on the conformation probability associated with
    the given complex_name. complex_name should match the name field of one of
    the complex in the resting set used to construct this NupackSampleJob.
    The Bayesian estimator for the standard error is used
      SE = (n+1)*(N-n+1)/((N+3)*(N+2)*(N+2))
    to avoid artifacts when n=0,N.
    """
    index = self.get_complex_index(complex_name)
    Nc = self._complex_counts[index]
    N = self.total_sims
    return math.sqrt((Nc+1.0)*(N-Nc+1)/((N+3)*(N+2)*(N+2)))

  def get_complex_prob_data(self, complex_name = None):
    """
    Returns raw sampling data for the given complex_name. Data is returned as a
    list consisting of float values, where each float value is 1 minus the
    maximum fractional defect for any domain in that sampled secondary
    structure.
    """
    return self._data[complex_name]

  def set_complex_prob_data(self, complex_name, data):
    """
    Set the raw sampling data for the given complex_name. Should be used only
    when importing an old KinDA session to restore state.
    """
    self._data[complex_name] = np.array(data)

  def get_num_sims(self):
    """
    Returns the total number of sampled secondary structures.
    """
    return self.total_sims

  def get_complex_count(self, complex_name = None):
    # TODO: fix me, why is the int necessary here if complex_name = None???
    return int(self._complex_counts[self.get_complex_index(complex_name)])

  def set_similarity_threshold(self, similarity_threshold):
    """
    Modifies the similarity threshold used to classify each sampled secondary
    structure as one of the expected complex conformations. The complex
    probability estimates are recomputed based on the new similarity threshold,
    so that future calls to get_complex_prob() will give correct statistics
    based on previously sampled secondary structures.
    """
    ## Change internal similarity threshold
    self.similarity_threshold = similarity_threshold

    ## Recalculate complex counts for new similarity threshold
    self.update_complex_counts()

  def sample(self, num_samples, status_func=None):
    """
    Calls sample_multiprocessing or sample_singleprocessing depending on the
    value of self.multiprocessing.
    """
    if self.multiprocessing:
      self.sample_multiprocessing(num_samples, status_func=status_func)
    else:
      self.sample_singleprocessing(num_samples, status_func=status_func)

  def sample_multiprocessing(self, num_samples, status_func=None):
    """
    Runs sample() in multiple processes.
    """
    # Temporarily remove SIGINT event handler from current process
    sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)

    # Create Pool instance and worker processes
    k = min(self._mp_ctx.cpu_count(), num_samples)
    p = self._mp_ctx.Pool(processes=k)

    # Restore original SIGINT event handler (if possible)
    if sigint_handler is None: sigint_handler = signal.SIG_DFL
    signal.signal(signal.SIGINT, sigint_handler)

    # Setup args for each process
    samples_per_worker = int(num_samples / k)
    args = [(self, samples_per_worker+1)] * (num_samples % k)
    args += [(self, samples_per_worker)] * (k - (num_samples % k))
    it = p.imap_unordered(sample_global, args)
    p.close()
    try:
      sims_completed = 0
      for cplx in it:
        self.add_sampled_complexes(cplx)
        sims_completed += len(cplx)
        if status_func is not None:
          status_func(sims_completed)
    except KeyboardInterrupt:
      print("\nSIGINT: Ending NUPACK sampling prematurely...")
      p.terminate()
      p.join()
      raise KeyboardInterrupt

  def sample_singleprocessing(self, num_samples, status_func=None):
    """
    Queries Nupack for num_samples secondary structures, sampled from the
    Boltzmann distribution of secondary structures for this resting set. The
    sampled secondary structures are automatically processed to compute new
    estimates for each conformation probability. The nupack_params dict that was
    given to this job during initialization is passed along to the Nupack Python
    interface.
    """
    results = sample_global((self, num_samples))
    self.add_sampled_complexes(results)
    if status_func is not None:
      status_func(len(results))

  def add_sampled_complexes(self, sampled):
    """
    Processes a list of sampled Complex objects, computing the similarity to
    each of the resting set conformations and updating the estimates for each
    conformation probability.
    """

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

  def update_complex_counts(self):
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
      
  def reduce_error_to(self, rel_goal, max_sims, complex_name = None, 
      init_batch_size = 50, 
      min_batch_size = 50, 
      max_batch_size = 1000,
      verbose = 0):
    """Stochastic sampling of secondary structures until the error-bars are satisified.

    Querys NUPACK for secondary structures sampled from the Boltzmann distribution
    until the standard error of the estimated probability for the given complex_name
    is at most rel_goal*complex_prob, or max_sims conformations have been sampled.
    If no complex_name is given, the halting condition is based on the error
    for the spurious conformation probability.

    Args:
      rel_goal (float): The realtive error goal.
      max_sims (int): Maximal number of sampled conformations.
      complex_name (str, optional): Standard error for sampling this complex.
        Defaults to None, the probability of a spurious conformation.
      init_batch_size (int, optional): Batch size for sampling when no data 
        has been collected.
      min_batch_size (int, optional): Minimum batch size for sampling before
        new error-bars are calculated.  
      max_batch_size (int, optional): Maximum batch size for sampling before 
        new error-bars are calculated.
      verbose (int, optional): Print a progress table. 0: silent mode,
        1: print the rows of a table. 2: print header and rows of a table,
        3: start a new row for every new batch. 4: start a new row whenever
        there is new data available. Defaults to 0.
    """
    def status_func(batch_sims_done, inline=True):
      # Update only the right part of the separator. We don't want to bias the
      # left side with temporary results from fast simulations.
      if verbose > 3: inline = False
      total_sims = self.total_sims
      total_success = self.get_complex_count(complex_name)
      total_failure = total_sims - total_success

      if exp_add_sims is None:
        assert num_sims + batch_sims_done == self.total_sims
        update_func([complex_name, prob, error, goal, " |", 
          "{:d}/{:d}".format(batch_sims_done,num_trials), 
          "{:d}/--".format(total_sims), 
          "{:d}/{:d}".format(total_success, total_failure), 
          "      {}".format('--')], inline) 
      else:
        update_func([complex_name, prob, error, goal, " |", 
            "{:d}/{:d}".format(batch_sims_done, num_trials), 
            "{:d}/{:d}".format(total_sims, max(0,exp_add_sims-batch_sims_done)), 
            "{:d}/{:d}".format(total_success, total_failure), 
            "est{:4d}%".format(100*total_sims/
              (total_sims+max(0,exp_add_sims-batch_sims_done)))], inline) 
 
    # Get initial values
    num_sims = 0 # The number of finished simulations.
    prob = self.get_complex_prob(complex_name)
    error = self.get_complex_prob_error(complex_name)
    goal = rel_goal * prob

    if verbose:
      if verbose > 1:
        if self.multiprocessing:
          print(f'#    [MULTIPROCESSING ON] (over {self._mp_ctx.cpu_count()} cores)')
        else:
          print('#    [MULTIPROCESSING OFF]')

      # Prepare progress table
      update_func = print_progress_table(
          ["complex", "prob", "error", "err goal", " |", "batch sims", "done/needed", "S/F", "progress"],
          col_widths = [11, 10, 10, 10, 4, 15, 15, 12, 9], 
          col_format_specs = ['{}', '{:.6f}', '{:1.4g}', '{:1.4g}', '{}', '{}', '{}', '{}', '{}'],
          skip_header = True if verbose == 1 else False)
      update_func([complex_name, prob, error, goal, " |", "--/--", "--/--", "--/--", "--"])

    # Run simulations
    while not error <= goal and num_sims < max_sims:
      # Estimate additional trials based on inverse square root relationship
      # between error and number of trials
      if self.total_sims == 0:
        num_trials = init_batch_size
        exp_add_sims = None
      else:
        exp_add_sims = max(0, int(self.total_sims * ((error/goal)**2 - 1) + 1))
        num_trials = max(
            min(max_batch_size, exp_add_sims, max_sims - num_sims, self.total_sims + 1), 
            min_batch_size)
        
      # Query Nupack
      if verbose:
        status_func(0) 
      self.sample(num_trials, status_func=status_func if verbose else None)
      if verbose:
        status_func(num_trials, inline=(verbose <= 2))

      # Update estimates and goal
      num_sims += num_trials
      prob = self.get_complex_prob(complex_name)
      error = self.get_complex_prob_error(complex_name)
      goal = rel_goal * prob

    exp_add_sims = max(0, int(self.total_sims * ((error/goal)**2 - 1) + 1))

    if verbose: # Return result
      tot_sims = self.total_sims
      total_success = self.get_complex_count(complex_name)
      total_failure = tot_sims - total_success
      update_func([complex_name, prob, error, goal, " |", 
        "{:d}/{:d}".format(num_sims, max_sims), 
        "{:d}/{:d}".format(tot_sims, exp_add_sims), 
        "{:d}/{:d}".format(total_success, total_failure), 
        "{:7d}%".format(100*tot_sims/max([0,tot_sims+exp_add_sims]))], inline=False) 
    return num_sims

  def get_top_MFE_structs(self, num) -> List[Tuple[str, float]]:
    """
    NOTE: actually, these are sampled suboptimal structures.
    """
    strands = next(iter(self.restingset.complexes)).strands
    strand_seqs = [strand.sequence for strand in strands]
    energy_gap = 0.1
    struct_list = []
    while len(struct_list) < num:
      # Call Multistrand's Nupack wrapper
      struct_list = nupack.subopt(strand_seqs, energy_gap, **self._nupack_params)
      energy_gap += 0.5
    return [(s.structure.dp(), s.energy) for s in struct_list]

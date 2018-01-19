from ..imports import multistrandhome, dnaobjectshome

import math

# Import Multistrand
import multistrand.objects as MSObjects
from multistrand.options import Options as MSOptions
from multistrand.system import SimSystem as MSSimSystem

from dnaobjects import utils, io_Multistrand, Macrostate, RestingSet, Complex

from .. import options
from datablock import Datablock

from sim_utils import print_progress_table

if options.multistrand_params['multithreading']:
  try:
    import multiprocessing
    import signal
  except ImportError:
    print "Could not import multiprocessing package. Try turning multiprocessing off in options.py"
    raise

# GLOBALS
TRAJECTORY_MODE = 128
TRANSITION_MODE = 256
FIRST_PASSAGE_MODE = 16
FIRST_STEP_MODE = 48

EXACT_MACROSTATE = 0
BOUND_MACROSTATE = 1
DISASSOC_MACROSTATE = 2
LOOSE_MACROSTATE = 3
COUNT_MACROSTATE = 4

# Custom statistical functions:
def rate_mean_func(datablock):
  """Computes the expected rate given that the sampled reaction times
  follow an exponential distribution with a mean of 1/r. In this case,
  the correct estimate for the rate is the harmonic mean of the r's.
  If no data has been collected, returns NaN."""
  data = datablock.get_data()
  if len(data) > 0:
    times = [1.0 / r for r in data if r != 0]
    return len(times) / sum(times)
  else:
    return float('nan')
def rate_error_func(datablock):
  """Estimates the standard error for the rate estimate from the
  standard error for the measured reaction times. Based on local linearity
  of r=1/t, the error in the rates has the same proportion of the estimated
  rate as the error in the times.
  If 1 or fewer data points are collected, returns float('inf')."""
  data = datablock.get_data()
  n = len(data)
  if n > 1:
    times = [1.0 / r for r in data if r != 0]
    time_mean = sum(times) / n
    time_std = math.sqrt(sum([(t - time_mean)**2 for t in times]) / (n - 1))
    time_error = time_std / math.sqrt(n)
    # Based on estimated local linearity of relationship between t and r=1/t
    return time_error / time_mean**2
  else:
    return float('inf')
def k1_mean_func(datablock):
  data = datablock.get_data()
  n = len(data)
  n_s = len([x for x in data if x != 0])
  if n_s >= 1:
    mean = sum(data) / n_s
    return mean * (n_s + 1.0)/(n + 2.0)
  else:
    return (n_s + 1.0)/(n + 2.0)
#def k1_error_func(datablock):
def bernoulli_mean_func(datablock):
  # Expectation of the probability based on Bayesian analysis
  data = datablock.get_data()
  n = len(data)
  n_s = sum(data)
  return (n_s + 1.0) / (n + 2.0)
def bernoulli_std_func(datablock):
  return math.sqrt(datablock.get_mean() * (1 - datablock.get_mean()))
def bernoulli_error_func(datablock):
  data = datablock.get_data()
  n = len(data)
#  if n > 0:
#    return datablock.get_std() / math.sqrt(n)
#  else:
#    return float('inf')
  # Expectation of the standard error based on Bayesian analysis
  n_s = sum(data)
  return math.sqrt((n_s+1.0)*(n-n_s+1)/((n+3)*(n+2)*(n+2)));
    


# Global function for performing a single simulation, used for multithreading
def run_sims_global(params):
  multijob = params[0]
  num_sims = params[1]

  ms_options = multijob.create_ms_options(num_sims)
  MSSimSystem(ms_options).start()
  return ms_options


# MultistrandJob class definition
class MultistrandJob(object):
  """Represents a simulation job to be sent to Multistrand. Allows the
  calculation of reaction rates/times and error bars on calculated data.
  Depending on the output mode, different statistics are calculated on
  the trajectory results.
  This is the parent class to the more useful job classes that compile
  specific information for each job mode type."""
  
  def __init__(self, start_state, stop_conditions, sim_mode):
    self.ms_params = self.setup_ms_params(start_state = start_state,
                                          stop_conditions = stop_conditions,
                                          mode = sim_mode)
    
    self.datablocks["overall_time"] = Datablock()
    self.datablocks["overall_rate"] = Datablock(mean_func = rate_mean_func,
                                                error_func = rate_error_func)

    self.total_sims = 0

                                                
  def setup_ms_params(self, *args, **kargs):

    ## Convert DNAObjects to Multistrand objects
    if all(map(lambda x: isinstance(x, RestingSet), kargs['start_state'])):
      resting_sets = kargs['start_state']
      start_complexes = []
      boltzmann = True
    elif all(map(lambda x: isinstance(x, Complex), kargs['start_state'])):
      resting_sets = []
      start_complexes = kargs['start_state']
      boltzmann = False
    else:
      assert False, "Starting state must be all complexes or all resting sets"
    stop_conditions = kargs['stop_conditions']
    complexes = start_complexes

    ms_data = io_Multistrand.to_Multistrand(
        complexes = complexes,
        resting_sets = resting_sets,
        macrostates = stop_conditions
    )
    domains_dict = dict(ms_data['domains'])
    strands_dict = dict(ms_data['strands'])
    complexes_dict = dict(ms_data['complexes'])
    resting_sets_dict = dict(ms_data['restingstates'])
    macrostates_dict = dict(ms_data['macrostates'])

    ## Set ms_params with all parameters needed to create an MS Options object on the fly
    if boltzmann:
      start_state = [resting_sets_dict[rs] for rs in resting_sets]
    else:
      start_state = [complexes_dict[c] for c in start_complexes]

    params = {
        'start_state':        start_state,
        'dangles':            options.multistrand_params['dangles'],
        'simulation_time':    options.multistrand_params['sim_time'],
        'parameter_type':     options.multistrand_params['param_type'],
        'substrate_type':     options.multistrand_params['substrate_type'],
        'rate_method':        options.multistrand_params['rate_method'],
        'simulation_mode':    kargs['mode'],
        'temperature':        options.multistrand_params['temp'],
        'boltzmann_sample':   boltzmann,
        'stop_conditions':    [macrostates_dict[m] for m in stop_conditions],
        'join_concentration': options.multistrand_params['join_concentration']
    }

    return params
    
  def get_statistic(self, reaction, stat = 'rate'):
    return self.datablocks[reaction + "_" + stat].get_mean()
  def get_statistic_error(self, reaction, stat = 'rate'):
    return self.datablocks[reaction + "_" + stat].get_error()
  def get_statistic_data(self, reaction = 'overall', stat = 'rate'):
    return self.datablocks[reaction + "_" + stat].get_data()
  def set_statistic_data(self, reaction, stat, data):
    block = self.datablocks[reaction + "_" + stat]
    block.clear_data()
    block.add_data(data)

  

  def create_ms_options(self, num_sims):
    """ Creates a fresh MS Options object using the arguments in self.ms_params. """
      
    o = MSOptions(
        start_state     = self.ms_params['start_state'],
        dangles         = self.ms_params['dangles'],
        simulation_time = self.ms_params['simulation_time'],
        parameter_type  = self.ms_params['parameter_type'],
        substrate_type  = self.ms_params['substrate_type'],
        rate_method     = self.ms_params['rate_method']
    )
    o.simulation_mode    = self.ms_params['simulation_mode']
    o.temperature        = self.ms_params['temperature']
    o.boltzmann_sample   = self.ms_params['boltzmann_sample']
    o.stop_conditions    = self.ms_params['stop_conditions']
    o.join_concentration = self.ms_params['join_concentration']
    o.DNA23Metropolis()

    o.num_simulations = num_sims
    
    return o

  def run_simulations(self, num_sims, sims_per_update = 1, status_func = lambda:None):
    ## Run simulations using multithreading if specified
    if options.multistrand_params['multithreading']:
      self.run_sims_multithreaded(num_sims, sims_per_update, status_func)
    else:
      self.run_sims_singlethreaded(num_sims, sims_per_update, status_func)

  def run_sims_multithreaded(self, num_sims, sims_per_update = 1, status_func = lambda:None):
    """
    Sets up a multiprocessing.Pool object to run simulations concurrently over multiple subprocesses.
    The Pool object has number of worker processes equal to the number of cpus, as given by
      multiprocessing.cpu_count()
    Due to a Python bug, SIGINT events (e.g. produced by Ctrl-C) do not cause the worker processes to terminate
    gracefully, causing the result-handling loop to hang indefinitely and forcing the main process to be halted
    externally. This is resolved by removing the SIGINT handler before creating the worker processes, so
    that all workers ignore SIGINT, allowing the main process to handle SIGINT and terminate them without hanging.
    The original SIGINT handler is restored immediately after creating the worker processes (otherwise the
    main process would ignore SIGINT as well). Note that this workaround produces a short time interval in which
    all SIGINT signals are ignored.
    """

    # Temporarily remove the SIGINT event handler
    sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)

    # Create Pool instance and worker processes
    k = multiprocessing.cpu_count()
    p = multiprocessing.Pool( processes = k )

    # Restore original SIGINT event handler
    if sigint_handler is None:  sigint_handler = signal.SIG_DFL
    signal.signal(signal.SIGINT, sigint_handler)

    it = p.imap_unordered(run_sims_global, [(self, 2)] * num_sims)
    
    try:
      sims_completed = 0
      for res in it:
        self.process_results(res)

        sims_completed += 1
        if sims_completed % sims_per_update == 0:
          status_func(sims_completed)
      p.close()
    except KeyboardInterrupt:
      # (More) gracefully handle SIGINT by terminating worker processes properly
      # and then allowing SIGINT to be handled normally
      print "SIGINT: Terminating Multistrand processes prematurely..."
      p.terminate()
      p.join()
      raise KeyboardInterrupt


  def run_sims_singlethreaded(self, num_sims, sims_per_update = 1, status_func = None):
    sims_completed = 0
    while sims_completed < num_sims:
      sims_to_run = min(sims_per_update, num_sims - sims_completed)

      results = run_sims_global((self, sims_to_run))
      self.process_results(results)

      sims_completed += sims_to_run

      status_func(sims_completed)

  def process_results(self, ms_options):
    results = ms_options.interface.results
    times = [r.time for r in results]
    rates = [1.0/t for t in times if t != 0]
    
    self.datablocks["overall_time"].add_data(times)
    self.datablocks["overall_rate"].add_data(rates)

    self.total_sims = self.total_sims + len(results)

      
  
  def reduce_error_to(self, rel_goal, max_sims, reaction = 'overall', stat = 'rate', max_batch_size = 500):
    """Runs simulations to reduce the error to rel_goal*mean or until max_sims is reached."""
    def status_func(batch_sims_done):
      table_update_func([block.get_mean(), block.get_error(), goal, "", "%d/%d"%(batch_sims_done,num_trials), "%d/%d"%(num_sims+batch_sims_done,num_sims+exp_add_sims), str(100*(num_sims+batch_sims_done)/(num_sims+exp_add_sims))+"%"])

    tag = reaction + "_" + stat
    block = self.datablocks[tag]

    num_sims = 0
    error = block.get_error()
    goal = rel_goal * block.get_mean()

    # Check if any simulations are necessary
    if error <= goal or num_sims >= max_sims:
      return

    if options.multistrand_params['multithreading']:
      print "[MULTIPROCESSING ON] (with %d cores)" % multiprocessing.cpu_count()
    else:
      print "[MULTIPROCESSING OFF]"

    table_update_func = print_progress_table(
        [stat, "error", "err goal", "", "batch sims", "overall sims", "progress"],
        [14, 14, 14, 6, 17, 17, 10])
    table_update_func([block.get_mean(), error, goal, "", "--/--", "--/--", "--"])
    while not error <= goal and num_sims < max_sims:
      # Estimate additional trials based on inverse square root relationship
      # between error and number of trials
      if error == float('inf'):
        num_trials = 5
        exp_add_sims = max_sims
      else:
        reduction = error / goal
        exp_add_sims = int(block.get_num_points() * (reduction**2 - 1) + 1)
        num_trials = min(exp_add_sims, max_batch_size, max_sims - num_sims, self.total_sims + 1)
        
      self.run_simulations(num_trials, sims_per_update = 10, status_func = status_func)

      num_sims += num_trials
      error = block.get_error()
      goal = rel_goal * block.get_mean()

    

class FirstPassageTimeModeJob(MultistrandJob):
  
  datablocks = {}
  tags = None
  
  ms_params = {}
  
  def __init__(self, start_state, stop_conditions):
      
    super(FirstPassageTimeModeJob, self).__init__(start_state,
                                                  stop_conditions,
                                                  FIRST_PASSAGE_MODE)
      
    self.tags = [sc.tag for sc in self.ms_params['stop_conditions']]
    for sc in self.ms_params['stop_conditions']:
      self.datablocks[sc.tag + "_time"] = Datablock()
      self.datablocks[sc.tag + "_rate"] = Datablock(mean_func = rate_mean_func,
                                                    error_func = rate_error_func)
  
  def process_results(self, ms_options):
    results = ms_options.interface.results
    
    times = [r.time for r in results]
    rates = [1.0/t for t in times if t != 0]
    self.datablocks["overall_time"].add_data(times)
    self.datablocks["overall_rate"].add_data(rates)
    
    for tag in self.tags:
      relevant_sims = filter(lambda x: x.tag == tag, results)
      times = [r.time for r in relevant_sims]
      rates = [1.0 / t for t in times if t != 0]
      self.datablocks[tag + "_time"].add_data(times)
      self.datablocks[tag + "_rate"].add_data(rates)

    self.total_sims = self.total_sims + len(results)
      
      
      
class TransitionModeJob(MultistrandJob):
  
  datablocks = {}
  states = None
  
  ms_params = {}
  
  def __init__(self, start_state, macrostates, stop_conditions):
    stop_conditions = macrostates[:]
    
    # This actually modifies the stop macrostates in place, so this could be bad
    for sc in stop_states:
      sc.name = "stop:" + sc.name
      stop_conditions.append(sc)
      
    super(TransitionModeJob, self).__init__(start_complex, stop_conditions, TRANSITION_MODE)
      
    self.states = [sc.tag for sc in self.ms_params['stop_conditions']]
  
  def get_statistic(self, start_states, end_states, stat = 'rate'):
    tag = self.get_tag(start_states, end_states)
    return self.datablocks[tag + "_" + stat].get_mean()
  def get_statistic_error(self, start_states, end_states, stat = 'rate'):
    tag = self.get_tag(start_states, end_states)
    return self.datablocks[tag + "_" + stat].get_error()
  
  def process_results(self, ms_options):
    results = ms_options.interface.results
    transition_paths = ms_options.interface.transition_lists
    
    times = [r.time for r in results]
    rates = [1.0/t for t in times if t != 0]
    self.datablocks["overall_time"].add_data(times)
    self.datablocks["overall_rate"].add_data(rates)
    
    new_data = {}
    for path in transition_paths:
      collapsed_path = collapse_transition_path(path)
      for start, end in zip(collapsed_path[0:-1], collapsed_path[1:]):
        time_diff = end[0] - start[0]
        key_start = ",".join(filter(lambda x: x[1], enumerate(start[1])))
        key_end = ",".join(filter(lambda x: x[1], enumerate(end[1])))
        key = key_start + "->" + key_end
        if key in new_data:
          new_data[key].append(time_diff)
        else:
          new_data[key] = [time_diff]
    
    for key, times in new_data.items():
      if (key + "_time") not in self.datablocks:
        self.datablocks[key + "_time"] = Datablock()
        self.datablocks[key + "_rate"] = Datablock(mean_func = rate_mean_func,
                                                   error_func = rate_error_func)
      self.datablocks[key + "_time"].add_data(times)
      self.datablocks[key + "_rate"].add_data([1.0/t for t in times if t != 0])

    self.total_sims = self.total_sims + len(results)
    
    
  def collapse_transition_path(transition_path):
    """transition path is a list of the form
       [[time1, [in_state1, in_state2, ...]]
        [time2, [in_state1, in_state2, ...]]
        ...]"""
    return filter(lambda x: sum(x[1])>0, transition_path)
    
  
  def reduce_error_to(self, rel_goal, max_sims, start_states, end_states, stat = 'rate'):
    super(TransitionModeJob, self).reduce_error_to(rel_goal, max_sims,
        self.get_tag(start_states, end_states), stat)
    
  def get_tag(self, start_states, end_states):
    assert all([s in self.states for s in start_states]), "Unknown start state given in %s" % start_states
    assert all([s in self.states for s in end_states]), "Unknown end state given in %s" % end_states
    key_start = ",".join([i for i,s
                            in enumerate(self.states)
                            if s in start_states])
    key_end = ",".join([i for i,s
                            in enumerate(self.states)
                            if s in end_states])
    return key_start + "->" + key_end


class FirstStepModeJob(MultistrandJob):
  
  datablocks = {}
  tags = None
  
  ms_params = {}
  
  def __init__(self, start_state, stop_conditions):
  
    super(FirstStepModeJob, self).__init__(start_state, stop_conditions, FIRST_STEP_MODE)
      
    self.tags = [sc.tag for sc in self.ms_params['stop_conditions']]
    self.tags.append("None")
    for tag in self.tags:
      self.datablocks[tag + "_prob"] = Datablock(mean_func = bernoulli_mean_func,
                                                  std_func = bernoulli_std_func,
                                                  error_func = bernoulli_error_func)
      self.datablocks[tag + "_kcoll"] = Datablock(mean_func = rate_mean_func,
                                                  error_func = rate_error_func)
      self.datablocks[tag + "_k1"] = Datablock(mean_func = k1_mean_func)
      self.datablocks[tag + "_k2"] = Datablock(mean_func = rate_mean_func,
                                                  error_func = rate_error_func)

  def process_results(self, ms_options):
    results = ms_options.interface.results
    
    times = [r.time for r in results]
    rates = [1.0/t for t in times if t != 0]
    self.datablocks["overall_time"].add_data(times)
    self.datablocks["overall_rate"].add_data(rates)
    
    for r in [r for r in results if r.tag == None]:
      r.tag = "None"

    for tag in self.tags:
      relevant_sims = filter(lambda x: x.tag == tag, results)
      successes = map(lambda x: int(x.tag == tag), results)
      kcolls = [r.collision_rate for r in relevant_sims]
      k1s = map(lambda x: int(x.tag == tag) * x.collision_rate, results)
      k2s = [1.0 / r.time for r in relevant_sims if r.time != 0]

      self.datablocks[tag + "_prob"].add_data(successes)
      self.datablocks[tag + "_kcoll"].add_data(kcolls)
      self.datablocks[tag + "_k1"].add_data(k1s)
      self.datablocks[tag + "_k2"].add_data(k2s)

    self.total_sims = self.total_sims + len(results)
      
    

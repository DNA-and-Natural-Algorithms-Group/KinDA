import math

# Import Multistrand
import multistrand.objects as MSObjects
from multistrand.options import Options as MSOptions
from multistrand.system import SimSystem as MSSimSystem

from ..objects import utils, io_Multistrand, Macrostate, RestingSet, Complex

from .. import options
#from datablock import Datablock

import sim_utils

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
    
    self._stats_funcs = {
        'time': (sim_utils.time_mean, sim_utils.time_std, sim_utils.time_error),
        'rate': (sim_utils.rate_mean, sim_utils.rate_std, sim_utils.rate_error)
    }

    self._ms_results = {'tags': [], 'times': []}

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
    return self._stats_funcs[stat][0](reaction, self._ms_results)
  def get_statistic_error(self, reaction, stat = 'rate'):
    return self._stats_funcs[stat][2](reaction, self._ms_results)

  def get_simulation_data(self):
    return self._ms_results
  def set_simulation_data(self, ms_results):
    self._ms_results = ms_results
  

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

    self._ms_results['tags'].extend([('overall',)*len(times)])
    self._ms_results['times'].extend(times)

    self.total_sims = self.total_sims + len(results)

  def reduce_error_to(self, rel_goal, max_sims, reaction = 'overall', stat = 'rate', max_batch_size = 500):
    """Runs simulations to reduce the error to rel_goal*mean or until max_sims is reached."""
    def status_func(batch_sims_done):
      table_update_func([calc_mean(), calc_error(), goal, "", "%d/%d"%(batch_sims_done,num_trials), "%d/%d"%(num_sims+batch_sims_done,num_sims+exp_add_sims), str(100*(num_sims+batch_sims_done)/(num_sims+exp_add_sims))+"%"])
    def calc_mean():
      return self._stats_funcs[stat][0](reaction, self._ms_results)
    def calc_error():  
      return self._stats_funcs[stat][2](reaction, self._ms_results)

#    tag = reaction + "_" + stat
#    block = self.datablocks[tag]

    num_sims = 0
#    error = block.get_error()
#    goal = rel_goal * block.get_mean()
    error = calc_error()
    goal = rel_goal * calc_mean()

    # Check if any simulations are necessary
    if error <= goal or num_sims >= max_sims:
      return

    if options.multistrand_params['multithreading']:
      print "[MULTIPROCESSING ON] (with %d cores)" % multiprocessing.cpu_count()
    else:
      print "[MULTIPROCESSING OFF]"

    table_update_func = sim_utils.print_progress_table(
        [stat, "error", "err goal", "", "batch sims", "overall sims", "progress"],
        [14, 14, 14, 6, 17, 17, 10])
#    table_update_func([block.get_mean(), error, goal, "", "--/--", "--/--", "--"])
    table_update_func([calc_mean(), error, goal, "", "--/--", "--/--", "--"])
    while not error <= goal and num_sims < max_sims:
      # Estimate additional trials based on inverse square root relationship
      # between error and number of trials
      if error == float('inf'):
        num_trials = 5
        exp_add_sims = max_sims
      else:
        reduction = error / goal
        exp_add_sims = int(self.total_sims * (reduction**2 - 1) + 1)
        num_trials = min(exp_add_sims, max_batch_size, max_sims - num_sims, self.total_sims + 1)
        
      self.run_simulations(num_trials, sims_per_update = 10, status_func = status_func)
      status_func(num_trials)
#      print self._ms_results

      num_sims += num_trials
#      error = block.get_error()
#      goal = rel_goal * block.get_mean()
      error = calc_error()
      goal = rel_goal * calc_mean()

    print

    

class FirstPassageTimeModeJob(MultistrandJob):
  
  def __init__(self, start_state, stop_conditions):
      
    super(FirstPassageTimeModeJob, self).__init__(start_state,
                                                  stop_conditions,
                                                  FIRST_PASSAGE_MODE)
      
    self.tags = [sc.tag for sc in self.ms_params['stop_conditions']]
#    for sc in self.ms_params['stop_conditions']:
#      self.datablocks[sc.tag + "_time"] = Datablock()
#      self.datablocks[sc.tag + "_rate"] = Datablock(mean_func = rate_mean_func,
#                                                    error_func = rate_error_func)
#  
  def process_results(self, ms_options):
    results = ms_options.interface.results
    
    tags = [('overall', x.tag) if x.tag!=None else ('None',) for x in results]
    times = [r.time for r in results]
    self._ms_results['tags'].extend(tags)
    self._ms_results['times'].extend(times)

    self.total_sims = self.total_sims + len(results)
      
      
class TransitionModeJob(MultistrandJob):
  
#  datablocks = {}
#  states = None
#  
#  ms_params = {}
#  
  def __init__(self, start_state, macrostates, stop_states):
    stop_conditions = macrostates[:]
    
    # This actually modifies the stop macrostates in place, so this could be bad
    for sc in stop_states:
      sc.name = "stop:" + sc.name
      stop_conditions.append(sc)
      
    super(TransitionModeJob, self).__init__(start_complex, stop_conditions, TRANSITION_MODE)
      
    self.states = [sc.tag for sc in self.ms_params['stop_conditions']]
  
  def get_statistic(self, start_states, end_states, stat = 'rate'):
    tag = self.get_tag(start_states, end_states)
    return self._stats_funcs[stat][0](tag, self._ms_results)
  def get_statistic_error(self, start_states, end_states, stat = 'rate'):
    tag = self.get_tag(start_states, end_states)
    return self._stats_funcs[stat][2](tag, self._ms_results)
  
  def process_results(self, ms_options):
    results = ms_options.interface.results
    transition_paths = ms_options.interface.transition_lists
    
    tags = [('overall', r.tag) if r.tag!=None else ('None',) for r in results]
    times = [r.time for r in results]
    self._ms_results['tags'].extend(tags)
    self._ms_results['times'].extend(times)
   
    new_data = {}
    for path in transition_paths:
      collapsed_path = collapse_transition_path(path)
      for start, end in zip(collapsed_path[0:-1], collapsed_path[1:]):
        time_diff = end[0] - start[0]
        key_start = ",".join(filter(lambda x: x[1], enumerate(start[1])))
        key_end = ",".join(filter(lambda x: x[1], enumerate(end[1])))
        key = key_start + "->" + key_end
        self._ms_results['tags'].append(key)
        self._ms_results['times'].append(time_diff)

    self.total_sims += len(results)
    
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
  
#  datablocks = {}
#  tags = None
#  
#  ms_params = {}
#  
  def __init__(self, start_state, stop_conditions):
  
    super(FirstStepModeJob, self).__init__(start_state, stop_conditions, FIRST_STEP_MODE)
      
    self.tags = [sc.tag for sc in self.ms_params['stop_conditions']]
    self.tags.append("None")
#    for tag in self.tags:
#      self.datablocks[tag + "_prob"] = Datablock(mean_func = bernoulli_mean_func,
#                                                  std_func = bernoulli_std_func,
#                                                  error_func = bernoulli_error_func)
#      self.datablocks[tag + "_kcoll"] = Datablock(mean_func = rate_mean_func,
#                                                  error_func = rate_error_func)
#      self.datablocks[tag + "_k1"] = Datablock(mean_func = k1_mean_func)
#      self.datablocks[tag + "_k2"] = Datablock(mean_func = rate_mean_func,
#                                                  error_func = rate_error_func)
    self._stats_funcs['prob'] = (sim_utils.bernoulli_mean, sim_utils.bernoulli_std, sim_utils.bernoulli_error)
    self._stats_funcs['kcoll'] = (sim_utils.kcoll_mean, sim_utils.kcoll_std, sim_utils.kcoll_error)
    self._stats_funcs['k1'] = (sim_utils.k1_mean, sim_utils.k1_std, sim_utils.k1_error)
    self._stats_funcs['k2'] = (sim_utils.k2_mean, sim_utils.k2_std, sim_utils.k2_error)

    self._ms_results['kcoll'] = []
    self._ms_results['k2'] = []

  def process_results(self, ms_options):
    results = ms_options.interface.results
    
    tags = [('overall', r.tag) if r.tag!=None else ('None',) for r in results]
    times = [r.time for r in results]
    kcolls = [r.collision_rate for r in results]
    self._ms_results['tags'].extend(tags)
    self._ms_results['times'].extend(times)
    self._ms_results['kcoll'].extend(kcolls)
    self._ms_results['k2'].extend([1.0/t if t!=0 else float('inf') for t in times])

    self.total_sims += len(results)
      
    

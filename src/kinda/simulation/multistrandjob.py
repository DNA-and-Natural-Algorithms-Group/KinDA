# multistrandjob.py
# Created by Joseph Berleant, 9-15-2014
#
# Defines classes for interfacing with Multistrand in its 4 main modes
# (trajectory, transition, first passage, and first step) and collecting
# and processing data relevant to each mode.

import itertools as it
import signal

import numpy as np
import multiprocess

from multistrand.options import Options as MSOptions
from multistrand.options import Literals as MSLiterals
from multistrand.system import SimSystem as MSSimSystem

from ..objects import io_Multistrand, RestingSet, Complex
from . import sim_utils


MS_TIMEOUT = MSLiterals.time_out
MS_NOINITIALMOVES = MSLiterals.no_initial_moves
MS_ERROR = MSLiterals.sim_error


def run_sims_global(job_spec):
  """Multiprocessing function for performing a single simulation.
  """
  (multijob, num_sims) = job_spec
  ms_options = multijob.create_ms_options(num_sims)
  MSSimSystem(ms_options).start()
  ms_options.free_sim_system()
  return ms_options

# MultistrandJob class definition
class MultistrandJob:
  """
  Represents a simulation job to be sent to Multistrand. Allows the calculation
  of reaction rates/times and error bars on calculated data. Depending on the
  output mode, different statistics are calculated on the trajectory results.
  This is the parent class to the more useful job classes that compile specific
  information for each job mode type.
  """
  def __init__(self, start_state, stop_conditions, sim_mode,
          boltzmann_selectors = None, 
          multiprocessing = True, 
          multistrand_params = {}):
    self._multistrand_params = dict(multistrand_params)
    self._ms_options_dict = self.setup_ms_params(
      start_state = start_state, stop_conditions = stop_conditions,
      mode = sim_mode, boltzmann_selectors = boltzmann_selectors)

    self.multiprocessing = multiprocessing
    self._mp_ctx = multiprocess.get_context('spawn')
    
    self._stats_funcs = {
        'time': (sim_utils.time_mean, sim_utils.time_std, sim_utils.time_error),
        'rate': (sim_utils.rate_mean, sim_utils.rate_std, sim_utils.rate_error)
    }

    self._tag_id_dict = {
      MS_TIMEOUT: -1,
      MS_NOINITIALMOVES: -2,
      MS_ERROR: -3,
      'overall': 0
    }
    self._ms_results = {
        'valid': np.array([], dtype=np.int8), 
        'tags': np.array([], np.int64), 
        'times': np.array([])
        }
    self._ms_results_buff = {
        'valid': np.array([], dtype=np.int8), 
        'tags': np.array([], np.int64), 
        'times': np.array([])
        }
    # extra information about invalid simulations, like timeouts
    self._ms_results_invalid = []

    self.total_sims = 0

  @property
  def multistrand_params(self):
    return dict(self._multistrand_params)

  @property
  def tag_id_dict(self):
    return self._tag_id_dict.copy()
                                                
  def setup_ms_params(self, *args, **kargs):
    ## Extract keyword arguments
    start_state = kargs['start_state']
    stop_conditions = kargs['stop_conditions']
    resting_sets = [x for x in start_state if isinstance(x, RestingSet)]
    complexes = [x for x in start_state if isinstance(x, Complex)]

    if kargs['boltzmann_selectors'] is None:
      boltzmann = False
      boltzmann_selectors = [None]*len(start_state)
    else:
      boltzmann = True
      boltzmann_selectors = kargs['boltzmann_selectors']

    ## Convert DNAObjects to Multistrand objects
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

    ms_start_state = [
            resting_sets_dict[elem]
            if (elem in resting_sets_dict)
            else complexes_dict[elem]
        for elem in kargs['start_state']
    ]
    ms_stop_conditions = list(it.chain(*[macrostates_dict[m] for m in stop_conditions]))

    ## Set boltzmann sampling status for all Multistrand start_state complexes
    for elem,boltzmann_func in zip(ms_start_state, boltzmann_selectors):
      elem.boltzmann_sample = boltzmann
      elem.sampleSelect = boltzmann_func

    ## Set ms_params with all parameters needed to create an MS Options object
    ## on the fly
    options_dict = dict(
        self._multistrand_params,
        start_state =         ms_start_state,
        simulation_mode =     kargs['mode'],
        stop_conditions =     ms_stop_conditions
    )

    return options_dict
    
  def get_statistic(self, reaction, stat = 'rate'):
    return self._stats_funcs[stat][0](self._tag_id_dict[reaction], self._ms_results)

  def get_statistic_error(self, reaction, stat = 'rate'):
    return self._stats_funcs[stat][2](self._tag_id_dict[reaction], self._ms_results)

  def get_simulation_data(self):
    return self._ms_results

  def set_simulation_data(self, ms_results):
    # copy data from ms_results to self._ms_results_buff, while preserving data
    # types of numpy arrays in self._ms_results_buff
    for k,v in ms_results.items():
      self._ms_results_buff[k].resize(len(v), refcheck=False)
      self._ms_results[k] = self._ms_results_buff[k]
      if len(v) > 0:
        np.copyto(self._ms_results_buff[k], v)
    self.total_sims = len(self._ms_results['tags'])

  def add_simulation_data(self, ms_results):
    # copy data from ms_results to self._ms_results_buff, while preserving data
    # types of numpy arrays in self._ms_results_buff
    dim = len(ms_results['tags'])

    for k,v in ms_results.items():
      self._ms_results_buff[k] = np.append(self._ms_results_buff[k], v)
      self._ms_results[k] = self._ms_results_buff[k]
    self.total_sims = len(self._ms_results['tags'])

  def get_invalid_simulation_data(self):
    return self._ms_results_invalid

  def set_invalid_simulation_data(self, invalid_sim_data):
    self._ms_results_invalid = invalid_sim_data[:]
  
  def create_ms_options(self, num_sims: int) -> MSOptions:
    """
    Creates a fresh MS Options object using the arguments in
    self._ms_options_dict.
    """
    return MSOptions(**dict(self._ms_options_dict, num_simulations = num_sims))

  def run_simulations(self, num_sims, sims_per_update=1, sims_per_worker=1,
                      status_func=None):
    ## Run simulations using multiprocessing if specified
    if self.multiprocessing:
      self.run_sims_multiprocessing(num_sims, sims_per_update, sims_per_worker, status_func)
    else:
      self.run_sims_singleprocessing(num_sims, sims_per_update, status_func)

  def run_sims_multiprocessing(self, num_sims, sims_per_update=1,
                               sims_per_worker=1, status_func=None):
    """
    Sets up a multiprocessing.Pool object to run simulations concurrently over
    multiple subprocesses.

    The Pool object has number of worker processes equal to the number of cpus,
    as given by `multiprocessing.cpu_count()`.`
    Due to a Python bug, SIGINT events (e.g. produced by Ctrl-C) do not cause
    the worker processes to terminate gracefully, causing the result-handling
    loop to hang indefinitely and forcing the main process to be halted
    externally. This is resolved by removing the SIGINT handler before creating
    the worker processes, so that all workers ignore SIGINT, allowing the main
    process to handle SIGINT and terminate them without hanging. The original
    SIGINT handler is restored immediately after creating the worker processes
    (otherwise the main process would ignore SIGINT as well). Note that this
    workaround produces a short time interval in which all SIGINT signals are
    ignored.
    """
    # Temporarily remove the SIGINT event handler
    sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)

    # Create Pool instance and worker processes
    k = min(self._mp_ctx.cpu_count(), num_sims)
    p = self._mp_ctx.Pool(processes=k)

    # Restore original SIGINT event handler
    if sigint_handler is None:  sigint_handler = signal.SIG_DFL
    signal.signal(signal.SIGINT, sigint_handler)

    # Setup args for each process
    args = [(self, sims_per_worker)] * (num_sims // sims_per_worker)
    if num_sims%sims_per_worker > 0: args += [(self, num_sims%sims_per_worker)]
    it = p.imap_unordered(run_sims_global, args)
    p.close()
    try:
      sims_completed = 0
      for res in it:
        self.process_results(res)
        sims_completed += len(res.interface.results)
        if status_func is not None and sims_completed % sims_per_update == 0:
          status_func(sims_completed)
    except KeyboardInterrupt:
      # (More) gracefully handle SIGINT by terminating worker processes properly
      # and then allowing SIGINT to be handled normally
      print("\nSIGINT: Terminating Multistrand processes prematurely...")
      p.terminate()
      p.join()
      raise KeyboardInterrupt

  def run_sims_singleprocessing(self, num_sims, sims_per_update=1, status_func=None):
    sims_completed = 0
    while sims_completed < num_sims:
      sims_to_run = min(sims_per_update, num_sims - sims_completed)

      results = run_sims_global((self, sims_to_run))
      self.process_results(results)

      sims_completed += sims_to_run
      if status_func is not None:
        status_func(sims_completed)

  def preallocate_batch(self, batch_size):
    self._ms_results_buff['valid'].resize(self.total_sims + batch_size, refcheck=False)
    self._ms_results_buff['tags'].resize(self.total_sims + batch_size, refcheck=False)
    self._ms_results_buff['times'].resize(self.total_sims + batch_size, refcheck=False)
    self._ms_results['valid'] = self._ms_results_buff['valid'][:self.total_sims]
    self._ms_results['tags'] = self._ms_results_buff['tags'][:self.total_sims]
    self._ms_results['times'] = self._ms_results_buff['times'][:self.total_sims]

  def process_results(self, ms_options):
    results = ms_options.interface.results
    n = len(results)

    times = [r.time for r in results]
    tags  = [self._tag_id_dict[r.tag] for r in results]
    valid = [t!=self._tag_id_dict[MS_TIMEOUT]
             and t!=self._tag_id_dict[MS_ERROR] for t in tags]

    self._ms_results_buff['valid'][self.total_sims:self.total_sims+n] = valid
    self._ms_results_buff['tags'][self.total_sims:self.total_sims+n] = tags
    self._ms_results_buff['times'][self.total_sims:self.total_sims+n] = times

    ## Store extra information about invalid simulations
    for i in range(n):
      if not valid[i]:
        self._ms_results_invalid.append({
            'simulation_index': i+self.total_sims,
            'end_time': times[i],
            'type': 'timeout' if tags[i]==self._tag_id_dict[MS_TIMEOUT] else 'error',
            'seed': results[i].seed,
            'start_structure': results[i].start_state,
            'end_state': [list(v) for v in ms_options.interface.end_states[i]]
        })

    self.total_sims = self.total_sims + n
    for k in self._ms_results:
      self._ms_results[k] = self._ms_results_buff[k][:self.total_sims]

  def reduce_error_to(self, rel_goal, max_sims,
      reaction = 'overall',
      stat = 'rate',
      init_batch_size = 100,
      min_batch_size = 50,
      max_batch_size = 1000,
      sims_per_update = 1,
      sims_per_worker = 1,
      verbose = 0):
    """Stochastic simulations to reduce the error to rel_goal*mean or until max_sims is reached.
    
    Args:
      rel_goal (float): The realtive error goal.
      max_sims (int): Maximal number of sampled conformations.
      reaction (string, optional): The tag of the reaction.
      stat (string, optional): The statistics to be collected. 
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
      # Show updated mean and error based on finished simulations. Do not update the 
      # Goal and exp_add_sims though, after all we don't want to bias these estimates
      # toward fast reactions at runtime ...
      if verbose > 3: inline = False
      total_sims = self.total_sims
      total_success = (self._ms_results['tags']==self._tag_id_dict[reaction]).sum()
      total_timeout = total_sims - int((self._ms_results['valid']).sum())
      total_failure = total_sims - total_success - total_timeout

      if exp_add_sims is None:
        table_update_func([mean, error, goal, " |",
          "{:d}/{:d}".format(batch_sims_done, num_trials),
          "{:d}/--".format(total_sims),
          "{}/{}/{}".format(total_success, total_failure, total_timeout),
          "{}".format('--')], inline)
      else :
        table_update_func([mean, error, goal, " |",
          "{:d}/{:d}".format(batch_sims_done, num_trials),
          "{:d}/{:d}".format(total_sims, max(0, exp_add_sims-batch_sims_done)),
          "{}/{}/{}".format(total_success, total_failure, total_timeout),
          "est {:.0%}".format(
            total_sims / (total_sims+max(0, exp_add_sims-batch_sims_done)))], inline)

    def calc_mean():
      return self._stats_funcs[stat][0](self._tag_id_dict[reaction], self._ms_results)

    def calc_error():
      return self._stats_funcs[stat][2](self._tag_id_dict[reaction], self._ms_results)

    num_sims = 0
    error = calc_error()
    mean = calc_mean()
    goal = rel_goal * mean

    if verbose:
      if verbose > 1:
        if self.multiprocessing:
          print(f'#    [MULTIPROCESSING ON] (over {self._mp_ctx.cpu_count()} cores)')
        else:
          print('#    [MULTIPROCESSING OFF]')

      val_fmt = "{:.4%}" if stat == "prob" else "{:.3e}"
      table_update_func = sim_utils.print_progress_table(
          [stat, "error", "err goal", " |", "batch sims", "done/needed", "S/F/T", "progress"],
          col_widths = [10, 10, 10, 4, 17, 17, 17, 9],
          col_format_specs = [val_fmt] * 3 + ['{}'] * 5)
      table_update_func([mean, error, goal, " |"
        "--/--", "--/--", "--/--/--", "--"])

    # Run simulations
    while (
        # await convergence criterion
        (error > goal and num_sims < max_sims)
        # force evaluation of first batch (e.g., even if `error-goal` is too high),
        # except when `max_sims == 0` (e.g., during `export_data()`)
        or (num_sims == 0 and max_sims > 0)):

      # Estimate additional trials based on inverse square root relationship
      # between error and number of trials
      if self.total_sims == 0 :
        num_trials = min(max_sims-num_sims, init_batch_size)
        exp_add_sims = None
      elif error == float('inf') or goal == 0.0:
        num_trials = min(max_sims-num_sims, max_batch_size)
        exp_add_sims = None
      else:
        exp_add_sims = int(self.total_sims * ((error / goal)**2 - 1) + 1)
        num_trials = max(min(max_batch_size, exp_add_sims, self.total_sims + 1), min_batch_size)
      num_trials = min(num_trials, max_sims - num_sims)
        
      self.preallocate_batch(num_trials)
      self.run_simulations(num_trials, 
          sims_per_update = sims_per_update, 
          sims_per_worker = sims_per_worker, 
          status_func = status_func if verbose else None)
      if verbose:
        status_func(num_trials, inline=(verbose <= 2))

      num_sims += num_trials
      error = calc_error()
      mean = calc_mean()
      goal = rel_goal * mean

    total_sims = self.total_sims
    total_success = (self._ms_results['tags']==self._tag_id_dict[reaction]).sum()
    total_timeout = total_sims - int((self._ms_results['valid']).sum())
    total_failure = total_sims - total_success - total_timeout

    if error == float('inf') or goal == 0.0:
      exp_add_sims = None
    else:
      exp_add_sims = max(0, int(total_sims * ((error/goal)**2 - 1) + 1))

    if verbose:
      if exp_add_sims is None:
        table_update_func([mean, error, goal, " |",
          "{:d}/{:d}".format(num_sims, max_sims),
          "{:d}/--".format(total_sims),
          "{:d}/{:d}/{:d}".format(total_success, total_failure, total_timeout),
          "{}".format('--')], inline=False)
      else :
        table_update_func([calc_mean(), error, goal, " |",
          "{:d}/{:d}".format(num_sims, max_sims),
          "{:d}/{:d}".format(total_sims, exp_add_sims),
          "{:d}/{:d}/{:d}".format(total_success, total_failure, total_timeout),
          "{:.1%}".format(total_sims / (total_sims+exp_add_sims))], inline=False)


class FirstPassageTimeModeJob(MultistrandJob):
  
  def __init__(self, start_state, stop_conditions, unimolecular_k1_scale = 1000, **kargs):
    super().__init__(start_state, stop_conditions, MSLiterals.first_passage_time,
                     **kargs)
    self.tags = [sc.tag for sc in self._ms_options_dict['stop_conditions']]
    self._tag_id_dict.pop('overall')
    self._tag_id_dict.update((t,i) for i,t in enumerate(sorted(self.tags)))
    
    self._stats_funcs['prob'] = (
      sim_utils.bernoulli_mean, sim_utils.bernoulli_std,
      sim_utils.bernoulli_error)
    self._stats_funcs['k1'] = (
        sim_utils.uni_k1_mean(unimolecular_k1_scale), # these are actually classes so they can be picklable
        sim_utils.uni_k1_std(unimolecular_k1_scale),  # I should probably rethink the whole organization
        sim_utils.uni_k1_error(unimolecular_k1_scale) # of all the stats functions
    )
    self._stats_funcs['k2'] = (
      sim_utils.uni_k2_mean, sim_utils.uni_k2_std, sim_utils.uni_k2_error)

  def process_results(self, ms_options):
    results = ms_options.interface.results
    n = len(results)
    
    tags = [self._tag_id_dict[r.tag] for r in results]
    valid = [t!=self._tag_id_dict[MS_TIMEOUT]
             and t!=self._tag_id_dict[MS_ERROR] for t in tags]
    times = [r.time for r in results]

    start_ind, end_ind = self.total_sims, self.total_sims+n
    self._ms_results_buff['valid'][start_ind:end_ind] = valid
    self._ms_results_buff['tags'][start_ind:end_ind] = tags
    self._ms_results_buff['times'][start_ind:end_ind] = times

    ## Store extra information about invalid simulations
    for i in range(n):
      if not valid[i]:
        self._ms_results_invalid.append({
            'simulation_index': i+self.total_sims,
            'end_time': times[i],
            'type': 'timeout' if tags[i]==self._tag_id_dict[MS_TIMEOUT] else 'error',
            'seed': results[i].seed,
            'start_structure': results[i].start_state,
            'end_state': [list(v) for v in ms_options.interface.end_states[i]]
        })

    self.total_sims = self.total_sims + n
    for k in self._ms_results:
      self._ms_results[k] = self._ms_results_buff[k][:self.total_sims]
      
      
class TransitionModeJob(MultistrandJob):
  ## Warning: this class is largely untested  
  def __init__(self, start_state, macrostates, stop_states, **kargs):
    stop_conditions = macrostates[:]
    
    # This actually modifies the stop macrostates in place, so this could be bad
    for sc in stop_states:
      sc.name = "stop:" + sc.name
      stop_conditions.append(sc)
      
    super().__init__(start_state, stop_conditions, MSLiterals.transition,
                     **kargs)
    self.states = [sc.tag for sc in self._ms_options_dict['stop_conditions']]
    self._tag_id_dict.pop('overall')
    self._tag_id_dict.update((t,i) for i,t in enumerate(sorted(set(self.states))))
  
  def get_statistic(self, start_states, end_states, stat = 'rate'):
    tag = self.get_tag(start_states, end_states)
    return self._stats_funcs[stat][0](self._tag_id_dict[tag], self._ms_results)

  def get_statistic_error(self, start_states, end_states, stat = 'rate'):
    tag = self.get_tag(start_states, end_states)
    return self._stats_funcs[stat][2](self._tag_id_dict[tag], self._ms_results)
  
  def process_results(self, ms_options):
    results = ms_options.interface.results
    n = len(results)
    transition_paths = ms_options.interface.transition_lists
    
    tags = [self._tag_id_dict[r.tag] for r in results]
    valid = [t!=self._tag_id_dict[MS_TIMEOUT]
             and t!=self._tag_id_dict[MS_ERROR] for t in tags]
    times = [r.time for r in results]

    ## Store extra information about invalid simulations
    for i in range(n):
      if not valid[i]:
        self._ms_results_invalid.append({
            'simulation_index': i+self.total_sims,
            'end_time': times[i],
            'type': 'timeout' if tags[i]==self._tag_id_dict[MS_TIMEOUT] else 'error',
            'seed': results[i].seed,
            'start_structure': results[i].start_state,
            'end_state': [list(v) for v in ms_options.interface.end_states[i]]
        })
 
    for path in transition_paths:
      collapsed_path = collapse_transition_path(path)
      for start, end in zip(collapsed_path[0:-1], collapsed_path[1:]):
        time_diff = end[0] - start[0]
        key_start = ",".join([x for x in enumerate(start[1]) if x[1]])
        key_end = ",".join([x for x in enumerate(end[1]) if x[1]])
        key = key_start + "->" + key_end
        if key not in self._tag_id_dict:
          self._tag_id_dict[key] = len(self._tag_id_dict)
        tags.append(self._tag_id_dict[key])
        times.append(time_diff)
        valid.append(True)

    # Make sure there's enough space
    self.preallocate_batch(self.total_sims+len(tags) - len(self._ms_results_buff['tags']))

    self._ms_results_buff['valid'][self.total_sims:self.total_sims+len(tags)] = valid
    self._ms_results_buff['tags'][self.total_sims:self.total_sims+len(tags)] = tags
    self._ms_results_buff['times'][self.total_sims:self.total_sims+len(times)] = times

    self.total_sims += len(results)
    for k in self._ms_results:
      self._ms_results[k] = self._ms_results_buff[k][:self.total_sims]
    
  def collapse_transition_path(transition_path):
    """transition path is a list of the form
       [[time1, [in_state1, in_state2, ...]]
        [time2, [in_state1, in_state2, ...]]
        ...]"""
    return [x for x in transition_path if sum(x[1])>0]
    
  
  def reduce_error_to(self, rel_goal, max_sims, start_states, end_states,
                      stat = 'rate', **kwargs):
    super().reduce_error_to(rel_goal, max_sims,
                            self.get_tag(start_states, end_states), stat, **kwargs)
    
  def get_tag(self, start_states, end_states):
    assert all([s in self.states for s in start_states]), "Unknown start state given in %s" % start_states
    assert all([s in self.states for s in end_states]), "Unknown end state given in %s" % end_states
    key_start = ",".join([i for i,s in enumerate(self.states) if s in start_states])
    key_end = ",".join([i for i,s in enumerate(self.states) if s in end_states])
    return key_start + "->" + key_end


class FirstStepModeJob(MultistrandJob):
  def __init__(self, start_state, stop_conditions, **kargs):
    super().__init__(start_state, stop_conditions, MSLiterals.first_step, **kargs)
    self._tag_id_dict.pop('overall')
    self._tag_id_dict.update(
      (t,i) for i,t in enumerate(sorted(set(
        sc.tag for sc in self._ms_options_dict['stop_conditions']))))
    self.tags = list(self._tag_id_dict.keys())

    self._stats_funcs['prob'] = (
      sim_utils.bernoulli_mean, sim_utils.bernoulli_std, sim_utils.bernoulli_error)
    self._stats_funcs['kcoll'] = (
      sim_utils.kcoll_mean, sim_utils.kcoll_std, sim_utils.kcoll_error)
    self._stats_funcs['k1'] = (
      sim_utils.k1_mean, sim_utils.k1_std, sim_utils.k1_error)
    self._stats_funcs['k2'] = (
      sim_utils.k2_mean, sim_utils.k2_std, sim_utils.k2_error)

    self._ms_results['kcoll'] = np.array([])
    self._ms_results_buff['kcoll'] = np.array([])

  def preallocate_batch(self, batch_size):
    self._ms_results_buff['valid'].resize(self.total_sims+batch_size, refcheck=False)
    self._ms_results_buff['tags'].resize(self.total_sims+batch_size, refcheck=False)
    self._ms_results_buff['times'].resize(self.total_sims+batch_size, refcheck=False)
    self._ms_results_buff['kcoll'].resize(self.total_sims+batch_size, refcheck=False)
    self._ms_results['valid'] = self._ms_results_buff['valid'][:self.total_sims]
    self._ms_results['tags'] = self._ms_results_buff['tags'][:self.total_sims]
    self._ms_results['times'] = self._ms_results_buff['times'][:self.total_sims]
    self._ms_results['kcoll'] = self._ms_results_buff['kcoll'][:self.total_sims]

  def process_results(self, ms_options):
    results = ms_options.interface.results
    n = len(results)

    tags = [self._tag_id_dict[r.tag] for r in results]
    valid = [t!=self._tag_id_dict[MS_TIMEOUT]
             and t!=self._tag_id_dict[MS_ERROR] for t in tags]
    times = [r.time for r in results]
    kcolls = [r.collision_rate for r in results]
    
    assert len(self._ms_results_buff['tags']) >= self.total_sims+n
    self._ms_results_buff['valid'][self.total_sims:self.total_sims+n] = valid
    self._ms_results_buff['tags'][self.total_sims:self.total_sims+n] = tags
    self._ms_results_buff['times'][self.total_sims:self.total_sims+n] = times
    self._ms_results_buff['kcoll'][self.total_sims:self.total_sims+n] = kcolls

    ## Store extra information about invalid simulations
    for i in range(n):
      if not valid[i]:
        self._ms_results_invalid.append({
            'simulation_index': i+self.total_sims,
            'end_time': times[i],
            'type': 'timeout' if tags[i]==self._tag_id_dict[MS_TIMEOUT] else 'error',
            'seed': results[i].seed,
            'start_structure': results[i].start_state,
            'end_state': [list(v) for v in ms_options.interface.end_states[i]]
        })

    self.total_sims += n
    for k in self._ms_results:
      self._ms_results[k] = self._ms_results_buff[k][:self.total_sims]

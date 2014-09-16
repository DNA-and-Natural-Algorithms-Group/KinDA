import math

# Import Multistrand
import multistrand.objects as MSObjects
from multistrand.options import Options as MSOptions
from multistrand.system import SimSystem as MSSimSystem

from datablock import Datablock

# GLOBALS
TRAJECTORY_MODE = 128
TRANSITION_MODE = 256
FIRST_PASSAGE_MODE = 16
FIRST_STEP_MODE = 48

# Custom statistical functions:
def rate_mean_func(datablock):
  data = datablock.get_data()
  times = [1.0 / t for t in data]
  time_mean = sum(times) / len(times)
  return 1.0 / time_mean
def rate_error_func(datablock):
  data = datablock.get_data()
  n = len(data)
  if n > 1:
    times = [1.0 / t for t in data]
    time_mean = sum(times) / n
    time_std = math.sqrt(sum([(t - time_mean)**2 for t in times]) / (n - 1))
    time_error = time_std / math.sqrt(n)
    # Based on estimated local linearity of relationship between t and r=1/t
    error = time_error / time_mean**2
  else:
    error = float('inf')
  return error
def bernoulli_std_func(datablock):
  return math.sqrt(datablock.get_mean() * (1 - datablock.get_mean()))
def bernoulli_error_func(datablock):
  data = datablock.get_data()
  n = len(data)
  return datablock.get_std() / math.sqrt(n)
    
# MultistrandJob class definition
class MultistrandJob(object):
  """Represents a simulation job to be sent to Multistrand. Allows the
  calculation of reaction rates/times and error bars on calculated data.
  Depending on the output mode, different statistics are calculated on
  the trajectory results."""
  datablocks = {}
  
  ms_options = None
  
  def __init__(self, ms_options):
    self.ms_options = ms_options
    
    self.datablocks["overall_time"] = Datablock()
    self.datablocks["overall_rate"] = Datablock(mean_func = rate_mean_func,
                                                error_func = rate_error_func)
    
  def get_overall_rate(self):
    return self.datablocks["overall_rate"].get_mean()
  def get_overall_rate_error(self):
    return self.datablocks["overall_rate"].get_error()
  def get_overall_time(self):
    return self.datablocks["overall_time"].get_mean()
  def get_overall_time_error(self):
    return self.datablocks["overall_time"].get_error()
  
  def run_simulations(self, num_sims):
    self.ms_options.num_simulations = num_sims
    print "Running %d simulations..." % num_sims
    MSSimSystem(self.ms_options).start()
    print "Processing %d simulations..." % num_sims
    results = self.ms_options.interface.results
    self.process_results()
    
  def process_results(self):
    results = self.ms_options.interface.results
    times = [r.time for r in results]
    rates = [1.0/t for t in times]
    
    self.datablocks["overall_time"].add_data(times)
    self.datablocks["overall_rate"].add_data(rates)
    
    del self.ms_options.interface.results[:]
      
  
  def reduce_rate_error_to(self, threshold, type, tag):
    """type is either 'absolute' or 'relative', where 'relative' indicates
    that the threshold indicates a decimal fraction of the estimated rate.
    threshold is a list of two elements, indicating the lower and upper errors
    acceptable."""
    def get_updated_error():
      block = self.datablocks[tag + "_rate"]
      if type == 'absolute': 
        error_goal = threshold
      else:
        error_goal = threshold * block.get_mean()
      return (block.get_error(), error_goal)
    cur_error, error_goal = get_updated_error()
    block = self.datablocks[tag + "_rate"]
      
    while cur_error > error_goal:
      # Estimate additional trials based on inverse square root relationship
      # between error and number of trials
      reduction = cur_error / error_goal
      num_trials = int(block.get_num_points() * (reduction**2 - 1) + 1)
      num_trials = min(num_trials, 50)
      self.run_simulations(num_trials)
      cur_error, error_goal = get_updated_error()
    
  def reduce_rate_error_to(self, threshold, type):
    """type is either 'absolute' or 'relative', where 'relative' indicates
    that the threshold indicates a decimal fraction of the estimated rate.
    threshold is a list of two elements, indicating the lower and upper errors
    acceptable."""
    self.reduce_rate_error_to(threshold, type, "overall")
    

class FirstPassageTimeModeJob(MultistrandJob):
  
  datablocks = {}
  tags = None
  
  ms_options = None
  
  def __init__(self, ms_options):
    assert (ms_options.simulation_mode == TRAJECTORY_MODE
      or ms_options.simulation_mode == FIRST_PASSAGE_MODE)
      
    super(FirstPassageTimeModeJob, self).__init__(ms_options)
      
    self.tags = [sc.tag for sc in ms_options.stop_conditions]
    for sc in ms_options.stop_conditions:
      self.datablocks[sc.tag + "_time"] = Datablock()
      self.datablocks[sc.tag + "_rate"] = Datablock(mean_func = rate_mean_func,
                                                    error_func = rate_error_func)
  
  def get_reaction_rate(self, tag):
    return self.datablocks[tag + "_rate"].get_mean()
  def get_reaction_rate_error(self, tag):
    return self.datablocks[tag + "_rate"].get_error()
  def get_reaction_time(self, tag):
    return self.datablocks[tag + "_time"].get_mean()
  def get_reaction_time_error(self, tag):
    return self.datablocks[tag + "_time"].get_error()
  
  def process_results(self):
    results = self.ms_options.interface.results
    
    times = [r.time for r in results]
    rates = [1.0/t for t in times]
    self.datablocks["overall_time"].add_data(times)
    self.datablocks["overall_rate"].add_data(rates)
    
    for tag in self.tags:
      relevant_sims = filter(lambda x: x.tag == tag, results)
      times = [r.time for r in relevant_sims]
      rates = [1.0 / t for t in times]
      self.datablocks[tag + "_time"].add_data(times)
      self.datablocks[tag + "_rate"].add_data(rates)
      
    del self.ms_options.interface.results[:]
      
      
class TransitionModeJob(MultistrandJob):
  
  datablocks = {}
  states = None
  
  ms_options = None
  
  def __init__(self, ms_options):
    assert ms_options.simulation_mode == TRANSITION_MODE
      
    super(TransitionModeJob, self).__init__(ms_options)
      
    self.states = [sc.tag for sc in ms_options.stop_conditions]
  
  def get_reaction_rate(self, start_states, end_states):
    tag = self.get_tag(start_states, end_states)
    return self.datablocks[tag + "_rate"].get_mean()
  def get_reaction_rate_error(self, start_states, end_states):
    tag = self.get_tag(start_states, end_states)
    return self.datablocks[tag + "_rate"].get_error()
  def get_reaction_time(self, start_states, end_states):
    tag = self.get_tag(start_states, end_states)
    return self.datablocks[tag + "_time"].get_mean()
  def get_reaction_time_error(self, start_states, end_states):
    tag = self.get_tag(start_states, end_states)
    return self.datablocks[tag + "_time"].get_error()
  
  def process_results(self):
    results = self.ms_options.interface.results
    transition_paths = self.ms_options.interface.transition_lists
    
    times = [r.time for r in results]
    rates = [1.0/t for t in times]
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
      self.datablocks[key + "_rate"].add_data([1.0/t for t in times])
    
    del self.ms_options.interface.results[:]
    del self.ms_options.interface.transition_lists[:]
    
  def collapse_transition_path(transition_path):
    """transition path is a list of the form
       [[time1, [in_state1, in_state2, ...]]
        [time2, [in_state1, in_state2, ...]]
        ...]"""
    return filter(lambda x: sum(x[1])>0, transition_path)
    
  
  def reduce_rate_error_to(self, threshold, type, start_states, end_states):
    """type is either 'absolute' or 'relative', where 'relative' indicates
    that the threshold indicates a decimal fraction of the calculated
    statistic."""
    self.reduce_rate_error_to(threshold, type, self.get_tag(start_states, end_states))
    
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
  
  ms_options = None
  
  def __init__(self, ms_options):
    assert ms_options.simulation_mode == FIRST_STEP_MODE
      
    super(FirstStepModeJob, self).__init__(ms_options)
      
    self.tags = [sc.tag for sc in ms_options.stop_conditions]
    self.tags.append("None")
    for tag in self.tags:
      self.datablocks[tag + "_prob"] = Datablock(std_func = bernoulli_std_func,
                                                  error_func = bernoulli_error_func)
      self.datablocks[tag + "_kcoll"] = Datablock(mean_func = rate_mean_func,
                                                  error_func = rate_error_func)
      self.datablocks[tag + "_k2"] = Datablock(mean_func = rate_mean_func,
                                                  error_func = rate_error_func)

  
  def get_reaction_prob(self, tag):
    return self.datablocks[tag + "_prob"].get_mean()
  def get_reaction_prob_error(self, tag):
    return self.datablocks[tag + "_prob"].get_error()
  def get_reaction_kcoll(self, tag):
    return self.datablocks[tag + "_kcoll"].get_mean()
  def get_reaction_kcoll_error(self, tag):
    return self.datablocks[tag + "_kcoll"].get_error()
  def get_reaction_k1(self, tag):
    return self.get_reaction_prob(tag) * self.get_reaction_kcoll(tag)
  def get_reaction_k1_error(self, tag):
    """Assumes success and kcoll are uncorrelated. This is not the case."""
    n = self.datablocks[tag + "_kcoll"].get_num_points()
    prob = self.get_reaction_prob(tag)
    prob_err = self.get_reaction_prob_error(tag)
    kcoll = self.get_reaction_kcoll(tag)
    kcoll_err = self.get_reaction_kcoll_error(tag)
    return math.sqrt(prob * kcoll_err +
                     prob_err * kcoll +
                     prob_err * kcoll_err) / math.sqrt(n)
  def get_reaction_k2(self, tag):
    return self.datablocks[tag + "_k2"].get_mean()
  def get_reaction_k2_error(self, tag):
    return self.datablocks[tag + "_k2"].get_error()
  
  def process_results(self):
    results = self.ms_options.interface.results
    
    times = [r.time for r in results]
    rates = [1.0/t for t in times]
    self.datablocks["overall_time"].add_data(times)
    self.datablocks["overall_rate"].add_data(rates)
    
    for tag in self.tags:
      relevant_sims = filter(lambda x: x.tag == tag, results)
      successes = map(lambda x: int(x.tag == tag), results)
      kcolls = [r.collision_rate for r in relevant_sims]
      k2s = [1.0 / r.time for r in relevant_sims]
      self.datablocks[tag + "_prob"].add_data(successes)
      self.datablocks[tag + "_kcoll"].add_data(kcolls)
      self.datablocks[tag + "_k2"].add_data(k2s)
      
    del self.ms_options.interface.results[:]
    
  def reduce_rate_error_to(self, threshold, type, tag):
    """type is either 'absolute' or 'relative', where 'relative' indicates
    that the threshold indicates a decimal fraction of the estimated rate.
    threshold is a list of two elements, indicating the lower and upper errors
    acceptable."""
    def get_updated_error():
      block = self.datablocks[tag]
      if type == 'absolute': 
        error_goal = threshold
      else:
        error_goal = threshold * block.get_mean()
      return (block.get_error(), error_goal)
    cur_error, error_goal = get_updated_error()
    block = self.datablocks[tag]
      
    while cur_error > error_goal:
      # Estimate additional trials based on inverse square root relationship
      # between error and number of trials
      print "Error should be reduced from %s to %s" % (cur_error, error_goal)
      reduction = cur_error / error_goal
      num_trials = int(block.get_num_points() * (reduction**2 - 1) + 1)
      num_trials = min(num_trials, 50)
      self.run_simulations(num_trials)
      cur_error, error_goal = get_updated_error()
      
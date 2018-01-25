# sim_utils.py
# Created by Joseph Berleant, 1/18/2018
#
# Misc utility functions used by simulation code

import sys
import math

def print_progress_table(col_headers, col_widths = None, col_init_data = None):
  """ Pretty prints a progress table to provide status updates when running simulations
  involving Nupack or Multistrand. Returns a progress update function that can be called
  with new column values to update the status shown to the user.
  col_headers is a list of header strings for each column in the table.
  col_widths (if given) is a list of ints giving the width of each column. The header
  strings and data strings are clipped to one less than the column width.
  col_init_data allows initial data to be displayed the first time the table is printed. 
  """

  def update_progress(col_data):
    print ' '.join([(str(d)+' '*(w-1))[:w-1] for d,w in zip(col_data, col_widths)]) + '\r',
    sys.stdout.flush()

  if col_widths is None:
    col_widths = [max(len(h)+1, 8) for h in col_headers]

  header = ' '.join([(h+' '*(w-1))[:w-1] for h,w in zip(col_headers, col_widths)])
  print header

  if col_init_data is not None:
    update_progress(col_init_data)

  return update_progress


################################
# CUSTOM STATISTICAL FUNCTIONS
################################ 

def time_mean(success_tag, ms_results):
  """ Returns the average success time of a simulation. """
  success_times = [time for tag,time in zip(ms_results['tags'],ms_results['times']) if success_tag in tag]
  if len(success_times) > 0:
    return sum(success_times) / len(success_times)
  else:
    return float('nan')
def time_std(success_tag, ms_results):
  """ Returns the sample standard deviation for a successful simulation time. """
  success_times = [time for tag,time in zip(ms_results['tags'],ms_results['times']) if success_tag in tag]
  n = len(success_times)
  if n > 1:
    mean = sum(success_times) / n
    return math.sqrt(sum([(t - mean)**2 for t in success_times]) / (n-1))
  else:
    return float('inf')
def time_error(success_tag, ms_results):
  """ Returns the standard error on the mean simulation time. """
  success_times = [time for tag,time in zip(ms_results['tags'],ms_results['times']) if success_tag in tag]
  n = len(success_times)
  if n > 1:
    mean = sum(success_times) / n
    std = math.sqrt(sum([(t - mean)**2 for t in success_times]) / (n-1))
    return std / math.sqrt(n)
  else:
    return float('inf')

def rate_mean(success_tag, ms_results):
  """Computes the expected rate given that the sampled reaction times
  follow an exponential distribution with a mean of 1/r. In this case,
  the correct estimate for the rate is the harmonic mean of the r's.
  If no data has been collected, returns NaN."""
  tags = ms_results['tags']
  times = ms_results['times']
  success_times = [time for tag,time in zip(tags,times) if success_tag in tag]
  if len(success_times) > 0:
    return len(success_times) / sum(success_times)
  else:
    return float('nan')
def rate_std(success_tag, ms_results):
  """ The standard deviation for the random variable k = 1/E[t] is the
  same as the standard error reported by rate_error() """
  return rate_error(success_tag, ms_results)
def rate_error(success_tag, ms_results):
  """Estimates the standard error for the rate estimate from the
  standard error for the measured reaction times. Based on local linearity
  of r=1/t, the error in the rates has the same proportion of the estimated
  rate as the error in the times.
  If 1 or fewer data points are collected, returns float('inf')."""
  tags = ms_results['tags']
  times = ms_results['times']
  success_times = [time for tag,time in zip(tags,times) if success_tag in tag]
  n = len(success_times)
  if n > 1:
    time_mean = sum(success_times) / n
    time_std = math.sqrt(sum([(t - time_mean)**2 for t in success_times]) / (n - 1))
    time_error = time_std / math.sqrt(n)
    # Based on estimated local linearity of relationship between t and r=1/t
    return time_error / time_mean**2
  else:
    return float('inf')

def kcoll_mean(success_tag, ms_results):
  """ Computes the expected kcoll rate, given the sampled kcolls
  from Multistrand trajectories.
  If no kcoll values have been collected for this reaction, returns NaN. """
  tags = ms_results['tags']
  kcolls = ms_results['kcoll']
  success_kcolls = [kcoll for tag,kcoll in zip(tags,kcolls) if success_tag in tag]
  n_s = len(success_kcolls)
  if n_s > 0:
    return sum(success_kcolls) / n_s
  else:
    return float('nan')
def kcoll_std(success_tag, ms_results):
  """ Computes the standard deviation on kcoll, given the sampled values.
  If less than 2 kcoll values have been collected for this reaction, returns float('inf'). """
  tags = ms_results['tags']
  kcolls = ms_rsults['kcoll']
  success_kcolls = [kcoll for tag,kcoll in zip(tags,kcolls) if success_tag in tag]
  n_s = len(success_kcolls)
  if n_s > 1:
    mean = sum(success_kcolls) / n_s
    return math.sqrt(sum([(kcoll - mean)**2 for kcoll in success_kcolls]) / (n-1))
  else:
    return float('inf')
def kcoll_error(success_tag, ms_results):
  """ Computes the standard error on the expected value of kcoll.  """
  tags = ms_results['tags']
  success_kcolls = [kcoll for tag,kcoll in zip(tags,ms_results['kcoll']) if success_tag in tag]
  n_s = len(success_kcolls)
  if n_s > 1:
    mean = sum(success_kcolls) / n_s
    return mean / math.sqrt(n_s - 1)
  else:
    return float('inf')

def k1_mean(success_tag, ms_results):
  """ Reports the expected value of k1, the rate constant for
  the bimolecular step of a resting-set reaction. """
#  print '***', success_tag
  tags = ms_results['tags']
  kcolls = ms_results['kcoll']
  success_kcolls = [kcoll for tag,kcoll in zip(tags,kcolls) if success_tag in tag]
  n = len(kcolls)
  n_s = len(success_kcolls)
  if n_s >= 1:
    mean = sum(success_kcolls) / n_s
    return mean * (n_s + 1.0)/(n + 2.0)
  else:
    return float('nan')
def k1_std(success_tag, ms_results):
  """ The standard deviation for k1 is the
  same as the standard error reported by k1_error() """
  return k1_error(success_tag, ms_results)
def k1_error(success_tag, ms_results):
  """ Reports the standard error on the expected value of k1.
  See the KinDA paper for a derivation. """
  tags = ms_results['tags']
  kcolls = ms_results['kcoll']
  success_kcolls = [kcoll for tag,kcoll in zip(tags,kcolls) if success_tag in tag]
  n = len(kcolls)
  n_s = len(success_kcolls)
  if n_s > 1:
    gamma = sum(success_kcolls)
    return gamma * math.sqrt((n_s+1.) / (n_s * (n+2.)) * ( (n_s+2.)/((n_s-1.)*(n+3.)) - (n_s+1.)/(n_s*(n+2.)) ))
  else:
    return float('inf')

def bernoulli_mean(success_tag, ms_results):
  """ Expectation of the bernoulli random variable S_i based on Bayesian analysis """
  tags = ms_results['tags']
  n = len(tags)
  n_s = len([t for t in tags if success_tag in t])
  return (n_s + 1.0) / (n + 2.0)
def bernoulli_std(success_tag, ms_results):
  """ Returns the standard deviation of the bernoulli random variable S_i
  which is 1 iff the reaction was successful
  (see KinDA paper for a more complete definition). """
  mean = bernoulli_mean(success_tag, ms_results)
  return math.sqrt(mean * (1 - mean))
def bernoulli_error(success_tag, ms_results):
  """ Returns the standard error of the mean of S_i. """
  tags = ms_results['tags']
  n = len(tags)
  n_s = len([t for t in tags if success_tag in t])
  return math.sqrt((n_s+1.0)*(n-n_s+1)/((n+3)*(n+2)*(n+2)));
 
def k2_mean(success_tag, ms_results):
  """ Returns the expectation of k2, the rate constant for the
  unimolecular step of a resting-set reaction. """
  tags = ms_results['tags']
  success_t2s = [1/k2 for tag,k2 in zip(tags, ms_results['k2']) if success_tag in tag]
  success_kcolls = [kcoll for tag,kcoll in zip(tags, ms_results['kcoll']) if success_tag in tag]
  if len(success_t2s) > 0:
    return sum(success_kcolls) / sum([kcoll*t2 for kcoll,t2 in zip(success_kcolls, success_t2s)])
  else:
    return float('nan')
def k2_std(success_tag, ms_results):
  """ The standard deviation is the same as the error on the estimate for k2. """
  return k2_error(success_tag, ms_results)
def k2_error(success_tag, ms_results):
  tags = ms_results['tags']
  success_t2s = [1/k2 for tag,k2 in zip(tags, ms_results['k2']) if success_tag in tag]
  if len(success_t2s) > 0:
    return math.sqrt(len(success_t2s)) / sum(success_t2s)
  else:
    return float('inf')

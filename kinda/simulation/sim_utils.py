# sim_utils.py
# Created by Joseph Berleant, 1/18/2018
#
# Misc utility functions used by simulation code

from __future__ import print_function

import sys
import math
import numpy as np

from multistrand.options import Options as MSOptions

def print_progress_table(col_headers, col_widths = None, col_init_data = None, col_format_specs = None):
  """ Live updates on progress with NUPACK and Multistrand computations.

  Note: This table has two rows. The 

  Args:
    col_headers (list(str)): The header of the table.
    col_widths (list(int), optional): Spacing of the table columns. Strings are
      clipped to width-1. (?)
    col_intit_data (list(), optional): Prints initial data into the first row.
    col_format_specs (list(), optional): ?

  Returns:
    A progress update function which overwrites the data row (or the last line on screen).
  """

  def update_progress(col_data, inline=True):
    """Print new data to your progress table."""
    str_data = [('{:<'+str(w-1)+'}').format(f.format(d))[:w-1] for d,w,f in zip(col_data, col_widths, col_format_specs)]
    print("#    {}{}".format(' '.join(str_data), "\r" if inline else "\n"), end='')
    sys.stdout.flush()

  if col_widths is None:
    col_widths = [max(len(h)+1, 8) for h in col_headers]
  else:
    assert len(col_widths) == len(col_headers)

  if col_format_specs is None:
    col_format_specs = ['{}'] * len(col_headers)
  else:
    assert len(col_format_specs) == len(col_headers)

  header = ' '.join([(h+' '*(w-1))[:w-1] for h,w in zip(col_headers, col_widths)])
  print("#    {}".format(header))

  if col_init_data is not None:
    update_progress(col_init_data)

  return update_progress


################################
# CUSTOM STATISTICAL FUNCTIONS
################################ 

def time_mean(success_tag, ms_results):
  """ Returns the average success time of a simulation. """
  success_times = np.ma.array(ms_results['times'], mask=(ms_results['tags']!=success_tag))
  n_s = np.sum(~success_times.mask)
  if n_s > 0:
    return float(success_times.mean())
  else:
    return float('nan')
def time_std(success_tag, ms_results):
  """ Returns the sample standard deviation for a successful simulation time. """
  success_times = np.ma.array(ms_results['times'], mask=(ms_results['tags']!=success_tag))
  n_s = np.sum(~success_times.mask)
  if n_s > 1:
    return float(success_times.std(ddof=1))
  else:
    return float('inf')
def time_error(success_tag, ms_results):
  """ Returns the standard error on the mean simulation time. """
  success_times = np.ma.array(ms_results['times'], mask=(ms_results['tags']!=success_tag))
  n_s = np.sum(~success_times.mask)
  if n_s > 1:
    return float(success_times.std(ddof=1) / math.sqrt(n_s))
  else:
    return float('inf')

def rate_mean(success_tag, ms_results):
  """Computes the expected rate given that the sampled reaction times
  follow an exponential distribution with a mean of 1/r. In this case,
  the correct estimate for the rate is the harmonic mean of the r's.
  If no data has been collected, returns NaN."""
  success_times = np.ma.array(ms_results['times'], mask=(ms_results['tags']!=success_tag))
  n_s = np.sum(~success_times.mask)
  if n_s > 0:
    return float(1./success_times.mean())
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
  success_times = np.ma.array(ms_results['times'], mask=(ms_results['tags']!=success_tag))
  n_s = np.sum(~success_times.mask)
  if n_s > 1:
    time_mean = success_times.mean()
    time_error = success_times.std(ddof=1) / math.sqrt(n_s)
    # Based on estimated local linearity of relationship between t and r=1/t
    return float(time_error / time_mean**2)
  else:
    return float('inf')

def kcoll_mean(success_tag, ms_results):
  """ Computes the expected kcoll rate, given the sampled kcolls
  from Multistrand trajectories.
  If no kcoll values have been collected for this reaction, returns NaN. """
  success_kcolls = np.ma.array(ms_results['kcoll'], mask=(ms_results['tags']!=success_tag))
  n = np.sum(ms_results['valid'])
  n_s = np.sum(~success_kcolls.mask)
  if n_s > 0:
    return float(success_kcolls.mean())
  elif n > 1:
    return float(ms_results['kcoll'].max())
  else:
    return float('nan')
def kcoll_std(success_tag, ms_results):
  """ Computes the standard deviation on kcoll, given the sampled values.
  If less than 2 kcoll values have been collected for this reaction, returns float('inf'). """
  success_kcolls = np.ma.array(ms_results['kcoll'], mask=(ms_results['tags']!=success_tag))
  n_s = np.sum(~success_kcolls.mask)
  if n_s > 1:
    return float(success_kcolls.std(ddof=1))
  else:
    return float('inf')
def kcoll_error(success_tag, ms_results):
  """ Computes the standard error on the expected value of kcoll.  """
  success_kcolls = np.ma.array(ms_results['kcoll'], mask=(ms_results['tags']!=success_tag))
  n_s = np.sum(~success_kcolls.mask)
  if n_s > 1:
    return float(success_kcolls.std(ddof=1) / math.sqrt(n_s))
  else:
    return float('inf')

def k1_mean(success_tag, ms_results):
  """ Reports the expected value of k1, the rate constant for
  the bimolecular step of a resting-set reaction. """
  success_kcolls = np.ma.array(ms_results['kcoll'], mask=(ms_results['tags']!=success_tag))
  n = np.sum(ms_results['valid'])
  n_s = np.sum(~success_kcolls.mask)
  tags = ms_results['tags']
  if n_s > 0:
    return float(np.sum(success_kcolls) / (n + 2.0))
  elif n > 0:
    return float(ms_results['kcoll'].max() * (n_s + 1.0) / (n + 2.0))
  else:
    return float('nan')
def k1_std(success_tag, ms_results):
  """ The standard deviation for k1 is the
  same as the standard error reported by k1_error() """
  return k1_error(success_tag, ms_results)
def k1_error(success_tag, ms_results):
  """ Reports the standard error on the expected value of k1.
  See the KinDA paper for a derivation. """
  success_kcolls = np.ma.array(ms_results['kcoll'], mask=(ms_results['tags']!=success_tag))
  n = np.sum(ms_results['valid'])
  n_s = np.sum(~success_kcolls.mask)
  if n_s > 0:
    gamma = np.sum(success_kcolls)
    return float(gamma/(n+2.) * math.sqrt((2.*n - n_s + 1.) / (n_s * (n+3.))))
  else:
    return float('inf')

def bernoulli_mean(success_tag, ms_results):
  """ Expectation of the bernoulli random variable S_i based on Bayesian analysis """
  tag_mask = ms_results['tags']==success_tag
  n = np.sum(ms_results['valid'])
  n_s = tag_mask.sum()
  return (n_s + 1.0) / (n + 2.0)
def bernoulli_std(success_tag, ms_results):
  """ Returns the standard deviation of the bernoulli random variable S_i
  which is 1 iff the reaction was successful
  (see KinDA paper for a more complete definition). """
  mean = bernoulli_mean(success_tag, ms_results)
  return math.sqrt(mean * (1 - mean))
def bernoulli_error(success_tag, ms_results):
  """ Returns the standard error of the mean of S_i. """
  tag_mask = ms_results['tags']==success_tag
  n = np.sum(ms_results['valid'])
  n_s = tag_mask.sum()
  return math.sqrt((n_s+1.0)*(n-n_s+1)/((n+3)*(n+2)*(n+2)));
 
def k2_mean(success_tag, ms_results):
  """ Returns the expectation of k2, the rate constant for the
  unimolecular step of a resting-set reaction. """
  success_kcolls = np.ma.array(ms_results['kcoll'], mask=(ms_results['tags']!=success_tag))
  success_t2s = np.ma.array(ms_results['times'], mask=(ms_results['tags']!=success_tag))
  n_s = np.sum(~success_t2s.mask)
  if n_s > 0:
    return float(np.sum(success_kcolls) / np.sum(success_kcolls*success_t2s))
  else:
    return float('nan')
def k2_std(success_tag, ms_results):
  """ The standard deviation is the same as the error on the estimate for k2. """
  return k2_error(success_tag, ms_results)
def k2_error(success_tag, ms_results):
  success_kcolls = np.ma.array(ms_results['kcoll'], mask=(ms_results['tags']!=success_tag))
  success_t2s = np.ma.array(ms_results['times'], mask=(ms_results['tags']!=success_tag))
  n_s = np.sum(~success_t2s.mask)
  if n_s > 1:
    time_mean = np.sum(success_kcolls*success_t2s) / np.sum(success_kcolls)
    time_std = math.sqrt(np.sum(success_kcolls * (success_t2s - time_mean)**2) / (np.sum(success_kcolls)-1))
    time_err = time_std / math.sqrt(n_s)
    return float(time_err / time_mean**2)
  else:
    return float('inf')

# sim_utils.py
# Created by Joseph Berleant, 1/18/2018
#
# Misc utility functions used by simulation code

import sys

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
    print '\r' + ' '.join([(str(d)+' '*(w-1))[:w-1] for d,w in zip(col_data, col_widths)]),
    sys.stdout.flush()

  if col_widths is None:
    col_widths = [max(len(h)+1, 8) for h in col_headers]

  header = ' '.join([(h+' '*(w-1))[:w-1] for h,w in zip(col_headers, col_widths)])
  print header

  if col_init_data is not None:
    update_progress(col_init_data)

  return update_progress

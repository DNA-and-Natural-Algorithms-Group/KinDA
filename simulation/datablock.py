import os, atexit
import random, math
import options
  
# Define default statistics functions    
def default_mean_func(datablock):
  """Default mean function"""
  data = datablock.get_data()
  if len(data) == 0:
    return float('nan')
  else:
    return sum(data) / len(data)  
def default_std_func(datablock):
  """Default standard deviation function"""
  data = datablock.get_data()
  n = len(data)
  mean = datablock.get_mean()
  return math.sqrt(sum([(val - mean)**2 for val in data]) / (n - 1))
def default_error_func(datablock):
  """Default error function"""
  data = datablock.get_data()
  n = len(data)
  if n > 1:
    mean = datablock.get_mean()
    std = datablock.get_std()
#    return [-std / math.sqrt(n) + mean, std / math.sqrt(n) + mean]
    return std / math.sqrt(n)
  else:
#    return [float("-inf"), float("inf")]
    return float('inf')
      

class Datablock(object):
  """Represents a set of data points sampled on an unknown
  quantity. Allows calculation of mean, standard deviation,
  and error bars, with optional functions to calculate each."""
  
  id_counter = 0

  id = None
  datatype = None
  
  data = []
  
  mean_func = None
  std_func = None
  error_func = None
  
  def __init__(self, datatype = 'float',
              mean_func = None, std_func = None, error_func = None):
    """Creates a new Datablock object with a unique
    id. Assigns the statistics functions, if given, or the default functions, if not."""
    
    # Assign unique id and increment global id counter
    self.id = Datablock.id_counter
    Datablock.id_counter += 1
      
    # Set datatype
    if datatype != 'float' and datatype != 'int' and datatype != 'string':
      print "WARNING: Invalid datatype %s" % datatype
      datatype = 'float'
    self.datatype = datatype
    
    # Assign default or given functions
    self.set_mean_func(mean_func)
    self.set_std_func(std_func)
    self.set_error_func(error_func)

  
  ## Data management
  def add_data(self, vals):
    """Appends the given values to data already stored with this block."""
    self.data.append(vals)
    
  def get_data(self):
    """Returns all data stored with this block. This may be a lot..."""
    return self.data
    
  def clear_data(self):
    self.data = []

  def export_data(self, path):
    f = open(path, 'w')
    print "Exporting data in datablock {0}".format(self.id)
    for d in self.data:
      f.write("{0}\n".format(d))
    f.close()
    print "Data export complete!"
  
  ## Data manipulation/statistics
  def get_num_points(self):
    return len(self.get_data())
  def get_mean(self):
    return self.mean_func(self)
  def get_std(self):
    return self.std_func(self)
  def get_error(self):
    return self.error_func(self)
  
  def set_mean_func(self, mean_func):
    """Sets the mean function to that provided, or
    sets it to the default if mean_func is None."""
    if mean_func == None:
      self.mean_func = default_mean_func
    else:
      self.mean_func = mean_func
  def set_std_func(self, std_func):
    """Sets the standard deviation function to that provided,
    or sets it to the default if std_func is None."""
    if std_func == None:
      self.std_func = default_std_func
    else:
      self.std_func = std_func
  def set_error_func(self, error_func):
    """Sets the error function to that provided, or
    sets it to the default if error_func is None."""
    if error_func == None:
      self.error_func = default_error_func
    else:
      self.error_func = error_func
  
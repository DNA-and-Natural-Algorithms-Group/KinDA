import os, atexit
import random, math
import options

blocklist = [] # Stores active blocks for deletion at program exit

# Define a cleanup function and register for execution when Python exits.
@atexit.register
def cleanup_all():
  for block in blocklist[:]:
    block.cleanup()

  
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
  
  id = None
  filepath = None
  datatype = None
  
  open = False
  data = None
  
  mean_func = None
  std_func = None
  error_func = None
  
  cleaned = False
  
  def __init__(self, datatype = 'float',
              mean_func = None, std_func = None, error_func = None):
    """Creates a new Datablock object with a unique
    id and creates the associated file. Assigns the error
    function, if given, or the default function, if none is given."""
    
    # Check that temp file directory exists
    base_dir = os.path.join(os.curdir, "temp")
    if not os.path.exists(base_dir):
      os.mkdir(base_dir)
    
    # Create unique unused id
    self.id = random.randrange(0,9999999)
    self.filepath = os.path.join(base_dir, "tmp" + str(self.id))
    while os.path.exists(self.filepath):
      self.id = random.randrange(0,9999999)
      self.filepath = os.path.join(base_dir, "tmp" + str(self.id))
      
    open(self.filepath, "w").close() # Create file
    
    # Set datatype
    if datatype != 'float' and datatype != 'int' and datatype != 'string':
      print "WARNING: Invalid datatype %s" % datatype
      datatype = 'float'
    self.datatype = datatype
    
    # Assign default or given functions
    self.set_mean_func(mean_func)
    self.set_std_func(std_func)
    self.set_error_func(error_func)
    
    # Add to blocklist (for deletion if Python exits badly)
    blocklist.append(self)
    
  
  def __del__(self):
    """Cleans up by calling cleanup()."""
    if not self.cleaned:
      self.cleanup()
  
  ## Opening/Closing
  #  Blocks may be opened and closed explicitly for improved efficiency.
  #  If a block is not explicitly opened, it will close automatically
  #  after each operation.
  #  Note that even in an open block, the data file will always be up-to-date.
  def open_block(self):
    if not self.open:
      self.data = self.read_data()
      self.open = True
  def close_block(self):
    self.data = None
    self.open = False
  def read_data(self):
    with open(self.filepath, "r") as f:
      if self.datatype == 'float':
        vals = [float(line) for line in f.readlines()]
      elif self.datatype == 'int':
        vals = [int(line) for line in f.readlines()]
      elif self.datatype == 'string':
        vals = [line.strip() for line in f.readlines()]
      else:
        assert False, "Invalid datatype %s" % self.datatype
    return vals
  
  ## Data management
  def add_data(self, vals):
    """Appends the given values to data already stored with this block."""
    with open(self.filepath, "a") as f:
      for v in vals:
        f.write(repr(v) + "\n")
    if self.open:
      data.append(vals)
    
  def get_data(self):
    """Returns all data stored with this block. This may be a lot..."""
    if self.open:
      return self.data
    else:
      return self.read_data()
    
  def clear_data(self):
    open(self.filepath, "w").close()
    if self.open:
      self.data = []
  
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
  
  def cleanup(self):
    """If options.flags['no-cleanup'] is not set, deletes
    the associated file. Otherwise, nothing is done.
    The datablock object should not be handled after calling this function"""
    if not options.flags['no-cleanup']:
      os.remove(self.filepath)
    blocklist.remove(self)
    self.cleaned = True
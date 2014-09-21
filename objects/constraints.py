from utils import base_group, base_complement
                  
def base_group_intersect(g1, g2):
  bases1 = set(base_group[g1])
  bases2 = set(base_group[g2])
  intersection = bases1.intersection(bases2)
  for key, val in base_group.items():
    if intersection == set(val):
      return key
  return "x"
  

class Constraints(str):
  """The Constraints class represents a set of constraints
  on a sequence of nucleotides, with common operations provided for
  finding the complement or intersection of constraints. A Constraints
  object is immutable."""
  def __init__(self, constraints):
    super(Constraints, self).__init__(constraints)
    
  @property
  def complement(self):
    c = "".join([base_complement[n] for n in reversed(self)])
    return Constraints(c)
    
  def intersection(self, other):
    if len(self) != len(other):
      print "Cannot intersect constraints {0} and {1}".format(self, other)
      return Constraints("")
      
    intersected = "".join([base_group_intersect(g1, g2)
                             for g1, g2
                             in zip(self, other)])
    return Constraints(intersected)

  def __add__(self, other):
    return self.intersection(other)
  
  
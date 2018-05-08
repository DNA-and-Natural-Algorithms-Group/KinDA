# Global DNA nucelotide groups
dna_base_group =  {"A": "A",   "T": "T",   "C": "C",   "G": "G",
                   "R": "AG",  "Y": "CT",  "W": "AT",  "S": "CG",  "M": "AC",  "K": "GT", 
                   "B": "CGT", "V": "ACG", "D": "AGT", "H": "ACT",
                   "N": "ACGT", "x": ""
                  }
dna_base_complement = {"A": "T",   "T": "A",   "C": "G",   "G": "C",
                       "R": "Y",   "Y": "R",   "W": "W",   "S": "S",   "M": "K",   "K": "M",
                       "B": "V",   "V": "B",   "D": "H",   "H": "D",
                       "N": "N",   "x": "x"
                      }
# Global RNA nucelotide groups
rna_base_group =  {"A": "A",   "U": "U",   "C": "C",   "G": "G",
                   "R": "AG",  "Y": "CU",  "W": "AU",  "S": "CG",  "M": "AC",  "K": "GU", 
                   "B": "CGU", "V": "ACG", "D": "AGU", "H": "ACU",
                   "N": "ACGU", "x": ""
                  }
rna_base_complement = {"A": "U",   "U": "A",   "C": "G",   "G": "C",
                       "R": "Y",   "Y": "R",   "W": "W",   "S": "S",   "M": "K",   "K": "M",
                       "B": "V",   "V": "B",   "D": "H",   "H": "D",
                       "N": "N",   "x": "x"
                      }
                  

# Check for material type (DNA vs. RNA) so we use the right nucleotide dictionaries
## NOTE: RNA is not implemented. Always uses DNA.
base_group = dna_base_group
base_complement = dna_base_complement



def base_group_intersect(g1, g2):
  bases1 = set(base_group[g1])
  bases2 = set(base_group[g2])
  intersection = bases1.intersection(bases2)
  for key, val in base_group.items():
    if intersection == set(val):
      return key
  return "x"
  

class Sequence(str):
  """The Sequence class represents a set of constraints
  on a sequence of nucleotides, with common operations provided for
  finding the complement or intersection of a constraints sequence. A Sequence
  object is immutable."""
  def __init__(self, sequence):
    super(Sequence, self).__init__(sequence)
    
  @property
  def complement(self):
    c = "".join([base_complement[n] for n in reversed(self)])
    return Sequence(c)
    
  def intersection(self, other):
    if len(self) != len(other):
      print "Cannot intersect constraint sequences {0} and {1}".format(self, other)
      return Sequence("")
      
    intersected = "".join([base_group_intersect(g1, g2)
                             for g1, g2
                             in zip(self, other)])
    return Sequence(intersected)

  def __add__(self, other):
    return self.intersection(other)
  
  

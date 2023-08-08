"""
structure.py

Copyright (c) 2014 Caltech. All rights reserved.
Coded by: Joseph Berleant (jberleant@dna.caltech.edu)

This module defines an abstract representation of a complex's
'structure' -- the binding relationship between all nucleotides
in the complex.

`structure`
    A specification for all bound base pairs in a complex.
"""

#####################
# PARSING UTILITIES #
#####################
    
def parse_dotparen(structure):
  """
  Parses the given dot-paren structure, returning a dict of the bound
  nucleotides.
  """
  strand_structs = structure.split("+")
  
  bond_dict = {}
  stack = []
  for strand_num, strand_struct in enumerate(strand_structs):
    i = 0
    for char in strand_struct:
      loc = (strand_num, i)
      if char == '.':  # nucleotide is unbound
        bond_dict[loc] = None
        i += 1
      elif char == '(':  # nucleotide is bound
        stack.append(loc)
        i += 1
      elif char == ')':  # nucleotide is bound to top element of stack
        j = stack.pop()
        bond_dict[loc] = j
        bond_dict[j] = loc
        i += 1
      elif char == '*' or char == '?':  # nucleotide binding is unspecified
        bond_dict[loc] = '?'
        i += 1
      elif char != ' ':
        raise ValueError("Invalid dot-paren structure.")

  return bond_dict

  
def parse_strandlist(structure):
  """
  Parses the given strand-list structure specification, returning a dict of the
  bonds in the structure.
  """
  bond_dict = {}
  for strand_num, strand_struct in enumerate(structure):
    for i, elem in enumerate(strand_struct):
      bond_dict[(strand_num, i)] = elem
  return bond_dict

    
def expand_domain_dotparen(structure, strands):
  """
  Expands a domain-level dot-paren specification into a nucleotide-level
  representation.
  """
  strand_structs = structure.replace(" ","").split("+")
  expanded = ""
  for strand, strand_struct in zip(strands, strand_structs):
    domains = strand.base_domains()
    for d, char in zip(domains, strand_struct):
      expanded += char * d.length
    expanded += "+"
  return expanded[:-1]

  
def expand_domain_strandlist(structure, strands):
  """
  Expands a domain-level strand-list specification into a nucleotide- level
  representation.
  """
  expanded = []
  for strand, strand_struct in zip(strands, structure):
    struct = []
    for d, bound in zip(strand.base_domains(), strand_struct):
      if bound == None:
        struct.extend([None] * d.length)
      else:
        b_strand_ind = bound[0]
        b_domain_ind = bound[1]
        b_domain_start = sum([bd.length
                                for bd
                                in strands[b_strand_ind].base_domains()[0:b_domain_ind]])
        b_domain_end = b_domain_start + d.length
        indices = reversed(list(range(b_domain_start, b_domain_end)))
        struct.extend(list(zip([bound[0]] * d.length, indices)))
    expanded.append(struct)
  return expanded
  
class Structure:
  """ A Structure object represents how the nucleotides in a complex
  are bound to each other. It may be input and output in a variety
  of forms, including dot-paren notation and strand-list notation
  (for Peppercorn). """
  def __init__(self, *args, **kargs):
    """
    Initialization:
    
    Keyword Arguments:
    structure [type=str, required]           -- Structure in dot-paren or
                                                strand-list notation.
    strands [type=list of strands, required] -- List of strands in this structure.
    """
    self._strands = kargs['strands']
    self._strand_order = list(range(len(self._strands)))
    
    self.structure = kargs['structure']

    self._dotparen = None
    self._strandlist = None
      
  @property
  def strands(self):
    """
    Returns the strands used in this structure, in the proper order.
    """
    return [self._strands[num] for num in self._strand_order]
  
  @property
  def structure(self):
    """
    Returns this structure in strand-list notation.
    """
    return self.to_strandlist()

  @structure.setter
  def structure(self, struct):
    """
    Sets this structure with a dot-paren string or strand-list representation.
    """
    self._dotparen = None # reset these as they are no longer valid
    self._strandlist = None

    domains_per_strand = [len(s.base_domains()) for s in self.strands]
    nucleos_per_strand = [s.length for s in self.strands]
    if isinstance(struct, list):  ## Assume strand-list notation
      # determine if domain- or nucleotide-level specification
      elems_per_strand = [len(s) for s in struct]
      if elems_per_strand == domains_per_strand:  ## Assume domain-level specification
        struct = expand_domain_strandlist(struct, self.strands)
      else:
        assert elems_per_strand == nucleos_per_strand, "Invalid strand-list structure."
      self._bond_dict = parse_strandlist(struct)
      self._pseudoknotted = False
    else:  ## Assume dot-paren notation
      # determine if domain- or nucleotide-level specification
      chars_per_strand = [elem.count("(") +
                          elem.count(")") +
                          elem.count(".") +
                          elem.count("*") for elem in struct.split("+")]
      if chars_per_strand == domains_per_strand:  ## Assume domain-level specification
        struct = expand_domain_dotparen(struct, self.strands)
      else:
        ## Check valid nucleotide-level specification
        assert chars_per_strand == nucleos_per_strand, "Invalid dot-paren structure."
      self._bond_dict = parse_dotparen(struct)
      self._pseudoknotted = self.check_pseudoknotted()
      
    # Since the structure now refers to the rotated strands, we need
    # to finalize this rotation.
    self._strands = self.strands
    self._strand_order = list(range(len(self._strands)))
    
  @property
  def pseudoknotted(self):
    """
    Returns if this structure is pseudoknotted under this strand ordering. The
    result is not really valid if the structure's binding is not completely
    specified (i.e. '?' appear in it)
    """
    return self._pseudoknotted
    
  ## Utility functions
  def bound_to(self, strand_num, index):
    """
    Returns the strand number and index of the nucleotide bound to the specified
    nucleotide, or None if it is unbound.
    """
    real_strand_num = self._strand_order[strand_num]
    real_bound_info = self._bond_dict[(real_strand_num, index)]
    if real_bound_info == None or real_bound_info == '?':
      return real_bound_info
    else:
      return (self._strand_order.index(real_bound_info[0]), real_bound_info[1])
    
  def to_dotparen(self):
    """
    Returns a representation of this structure in dot-paren notation. If the
    structure does not have a valid dot-paren representation (e.g. because of
    pseudoknots), it will fail.
    """
    assert not self._pseudoknotted, "Cannot express pseudoknotted structure with dot-paren."
    if self._dotparen is not None:  return self._dotparen

    dotparen = ""
    for strand_num, strand in enumerate(self.strands):
      for i in range(strand.length):
        bound = self.bound_to(strand_num, i)
        if bound == None:
          dotparen += '.'
        elif bound == '?':
          dotparen += '*'
        elif (strand_num, i) < bound:
          dotparen += '('
        else:
          dotparen += ')'
      dotparen += '+'

    self._dotparen = dotparen[:-1]
    return self._dotparen
    
  def to_strandlist(self):
    """
    Returns a representation of this structure in strand-list notation.
    """
    if self._strandlist is not None:  return self._strandlist

    strandlist = []
    for strand_num, strand in enumerate(self.strands):
      strand_struct = []
      for i in range(strand.length):
        strand_struct.append(self.bound_to(strand_num, i))
      strandlist.append(strand_struct)

    self._strandlist = strandlist
    return strandlist
    
  def check_pseudoknotted(self):
    """
    Checks if the structure may have pseudoknots (when strands are in their
    current order). This is not a technical pseudoknot, which would require
    checking other strand orders, but is suitable for most purposes.

    There may be a more efficient way to do this...
    """
    stack = []
    for strand_num, strand in enumerate(self.strands):
      for i in range(strand.length):
        bound = self.bound_to(strand_num, i)
        if bound != None and bound != '?' and (strand_num, i) < bound:
          stack.append((strand_num, i))
        elif bound != None and bound != '?' and bound < (strand_num, i):
          if bound != stack.pop():
            # There's a pseudoknot!
            return True
    return False
    
  ## Modifiers
  def rotate_strands(self, amount = 1):
    """
    Rotates the strands in this structure. Future queries into this structure
    will use the new strand ordering.
    """
    amount = amount % len(self._strands)
    self._strand_order = self._strand_order[amount:] + self._strand_order[:amount]
    self._dotparen = None
    self._strandlist = None

  def __str__(self):
    return self.to_dotparen()

  def __repr__(self):
    return self.to_dotparen()
    
  def __eq__(self, other):
    return self.to_dotparen() == other.to_dotparen()

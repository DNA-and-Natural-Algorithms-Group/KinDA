"""
complex.py

Copyright (c) 2010-2014 Caltech. All rights reserved.
Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
          Joseph Berleant (jberleant@dna.caltech.edu)

This module defines a simple dna `complex` object.

`complex`
    A completely connected set of nucleic acid strands.
"""


from strand import Strand
from structure import Structure

class Complex(object):
  """
  A representation of a single connected complex of strands.
  """
  id_counter = 0
  
  def __init__(self, *args, **kargs):
    """
    Initialization:
    
    Keyword Arguments:
    name [type=str]                -- Name of the complex. An automatic name
                                      "complex_###" is assigned if this is not
                                      given.
    strands [type=list of strands] -- List of strands in this complex
    structure [type=str OR list]   -- Structure for this complex, in dot-paren
                                      notation or strand-list notation
                                      as used by Casey's enumerator. Note that
                                      only the latter representation can be
                                      used for pseudoknotted structures.
    """
    # Assign id
    self.id = Complex.id_counter
    Complex.id_counter += 1
    
    # Assign DNA object type
    self._object_type = 'complex'
    
    # Assign name
    if 'name' in kargs: self.name = kargs['name']
    else: self.name = "complex_{0}".format(self.id)
    
    # Assign strands
    if 'strands' in kargs: self._strands = list(kargs['strands'])
    else: raise ValueError("'strands' key argument required for complex.")
    
    # Calculate length (the total number of nucleotides in the complex)
    self._length = sum([s.length for s in self.strands])
    
    # Assign structure
    if 'structure' in kargs: self.structure = kargs['structure']
    else: self.structure = '+'.join(['*'*s.length for s in self.strands])
  
  
  ## Basic properties
  @property
  def length(self):
    """ Returns the total number of nucleotides in the complex. """
    return self._length
    
  @property
  def structure(self):
    """ Returns the binding of this complex in a Structure object."""
    return self._structure
  @structure.setter
  def structure(self, new_struct):
    """ Sets the binding of this complex. Accepts a structure in dot-paren
    or strand-list notation. """
    self._structure = Structure(structure = new_struct,
                                strands = self._strands)
  def bound_to(self, strand_num, index):
    """ Returns the (strand_num, index) pair of the nucleotide bound to
    the specified nucleotide, or None if the specified nucleotide
    is unbound."""
    return self._structure.bound_to(strand_num, index)
    
  @property
  def pseduoknotted(self):
    """ Returns True if this complex's structure has a valid dot-paren
    representation. This corresponds to pseudoknottedness with this
    strand-ordering if there are no unspecified bonds in the structure."""
    return self._structure.pseudoknotted
    
  ## Equivalence
  def equivalent_to(self, other):
    """ Returns True if the two complexes are equivalent. Two complexes
    are equivalent if there is a strand rotation under which the lists of
    strands are equivalent and the structure is the same."""
    ### This is important but non-trivial to code...
    raise NotImplementedError("Cannot compare complexes for equivalence yet...")
    
    
  ## DNA object hierarchy
  @property
  def strands(self):
    """ Returns the ordered list of strands that compose this complex. """
    return self._strands
  def base_domains(self):
    """ Returns a list of the non-composite domain breakdown for each
    strand in the complex. """
    return [s.base_domains() for s in self._strands]
  
  
  ## (In)equality
  def __eq__(self, other): 
    """ Returns True if the two complexes have the same name. """
    return type(self)==type(other) and self.name == other.name
  def __ne__(self, other):
    """ Returns True if the complexes are not equal. """
    return not self.__eq__(other)
  def __hash__(self):
    """ Returns a hash value for this Complex. """
    return sum([ord(c) * pow(2,i)
      for i, c
      in enumerate(self._object_type + self.name)])
  
  
  ## Output
  def __str__(self):
    """ Human-readable output formatting for this Complex object. """
    strand_info = "[" + ", ".join([str(s) for s in self._strands]) + "]"
    return "Complex {0} {1}: {2}".format(self.name, strand_info, self.structure)
  def __repr__(self):
    return "Complex({0})".format(self.name)
    
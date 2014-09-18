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
                                      notation or strand-sequence list notation
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
    if 'strands' in kargs: self._strands = kargs['strands'][:]
    else: raise ValueError("'strands' key argument required for complex.")
    
    # Calculate length (this may not ever be useful...)
    self._length = sum([s.length for s in self.strands])
    
    # Assign structure
    if 'structure' in kargs: self.structure = kargs['structure']
    else: self.structure = ''
  
  ## Basic properties
  @property
  def length(self):
    return self._length
    
  @property
  def structure(self):
    ### TODO
  @structure.setter
  def structure(self, new_struct):
    ### TODO
  def bound_to(self, nucleotide_index):
    ### TODO
    
  @property
  def pseduoknotted(self):
    return self._pseudoknotted
    
  ## Equivalence
  def equivalent_to(self, other):
    ### TODO
    
    
  ## DNA object hierarchy
  @property
  def strands(self):
    return self._strands
  def base_domains(self):
    return [s.base_domains for s in self._strands]
  
  
  ## Output
  def __str__(self):
    ### TODO
    
    
  ## Auxiliary functions
  # Structure parsing
  def parse_dotparen(self, struct):
    ### TODO
  def parse_strandlist(self, struct):
    ### TODO
  
  def __str__( self ):
    return "\
Complex: {fieldnames[0]:>9}: '{0.name}'\n\
       : {fieldnames[1]:>9}: {0.sequence}\n\
       : {fieldnames[2]:>9}: {0.structure}\n\
       : {fieldnames[3]:>9}: {1}".format(
       self,
       [i.name for i in self.strand_list],
       fieldnames =('Name','Sequence','Structure','Strands') )

  def _init_parse_structure( self, structure ):
    strand_count = len(self.strand_list)
    base_count = sum(len(s.sequence) for s in self.strand_list)
    total_flat_length = base_count + strand_count - 1
   
    if len(structure) == total_flat_length:
      self._fixed_structure = structure
    else:
      domain_count = sum(len(s.domain_list) for s in self.strand_list)
      if len(structure) != domain_count + strand_count - 1:
        error_msg = "ERROR: Could not interpret the passed structure [{0}];".format(structure)
        if domain_count > 0:
          error_msg += "\
Expected a structure composed of characters from '.()+'\
and with either length [{0}] for a complete structure, or length [{1}]\
for a domain-level structure. If giving a domain-level structure, it \
should have the layout [{2}].".format( total_flat_length, 
              domain_count + strand_count - 1,
              "+".join(''.join('x'*len(d.sequence) for d in s.domain_list)\
                       for s in self.strand_list))
        else:
          error_msg += " Expected a complete structure with [{1}] total characters chosen from '.()+'.".format( total_flat_length )
        raise ValueError( error_msg )
      else:
        matched_list = zip( structure,
                            reduce( lambda x,y: x+[len(d.sequence) for d in y.domain_list] + [1], self.strand_list,[]))
                            # the reduce just composes the
                            # domain_lists into one big ordered list
                            # of domains
        self._fixed_structure = "".join(i[0]*i[1] for i in matched_list)

  def canonical_strand( self ):
    """Return the name of the `canonical` strand for this complex.

    This is either the strand name which appears first alphabetically,
    or if no strands are named, this is just the name of the complex.

    Return Value:
      -- The string containing the canonical name.
    """
    return min(self.strand_list, key=lambda x: x.name).name
  
  # def __len__(self):
  #   """ Length of a complex is the number of strands contained.

  #   Use the attribute :attr:`sequence_length` if you need the sequence length. """
  #   return len(self.strand_list)

  # @property
  # def sequence_length( self ):
  #   """ The total length of all contained strands. """
  #   return len("".join([i.sequence for i in self.strand_list]))


  @property
  def structure(self):
    """ If this complex is set to use boltzmann sampling, this property returns a newly sampled structure. Otherwise it gets the fixed structure for the complex."""
    return self._fixed_structure

  @structure.setter
  def structure(self,value):
    # I include the following due to the 'normal' use cases of Complex
    # involve it immediately being deepcopy'd (e.g. in starting
    # states, etc). So allowing someone to set the fixed structure
    # member is probably not a good idea - we could remove this
    # completely so the structure component is pretty much immutable.
    #
    # In practice, I can see a case where it may be useful to modify it if you're using an interactive mode, but even in those cases it may not do what the user wants, if they're relying on that changing previous parts. For example:
    # c1 = Complex("c1","c1", "......((.......)).....")
    # o = Options()           
    # o.start_state = [c1]
    # ... user runs a simulation and gets a syntax warning about mismatched parens...
    # c1.structure = "......((........))...."
    #
    # # This did NOT actually change o at all! So the warning used
    # # here mentions that, and just in case, returns a new object
    # # anyways!
    #
    import warnings
    warnings.warn("Setting a Complex's structure does not [usually] change existing uses of this Complex, so the object returned is a NEW object to avoid any confusion as to how it may affect previous usages.",SyntaxWarning)
    import copy
    retval = copy.deepcopy(self)
    retval._fixed_structure = value
    return retval
    
  @property
  def fixed_structure(self):
    """ The structure used to create this Complex. """
    return self._fixed_structure

  @property
  def sequence(self):
    """ The calculated 'flat' sequence for this complex. """
    return "+".join([strand.sequence for strand in self.strand_list])
  

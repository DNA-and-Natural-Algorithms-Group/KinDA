"""
strand.py

Copyright (c) 2010-2014 Caltech. All rights reserved.
Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
          Joseph Berleant (jberleant@dna.caltech.edu)

This module defines a simple dna `strand` object.

`strand`
    A nucleic acid strand, logically composed of either `domain` objects or nucleic acid bases; physically composed of a connected sequence of nucleic acid bases.
"""


from domain import Domain

class Strand(object):
  """
  Represents a Strand object, composed of a series of domains or a set of
  sequence constraints.
  """
  id_counter = 0

  def __init__(self, *args, **kargs ):
    """
    Initializes a new Strand object.  The following keyword arguments are
    accepted: name, domains OR constraints.
    
    If direct sequence constraints are specified, a new domain is defined
    with these constraints, and is assigned as the domain list for this
    Strand.
    """
    
    # Assign id
    self.id = Strand.id_counter
    Strand.id_counter += 1
    
    # Assign DNA object type
    self._object_type = 'strand'
    
    # Assign name
    if 'name' in kargs: self.name = kargs['name']
    else: self.name = 'strand_{0}'.format(self.id)
    
    # If sequence constraints were specified, create a dummy domain with
    # these constraints. Otherwise, assign the given list of domains.
    if 'domains' in kargs and 'constraints' not in kargs:
      self._domains = kargs['domains']
    elif 'constraints' in kargs and 'domains' not in kargs:
      d = Domain(name = self.name + "_domain",
                 constraints = kargs['constraints'])
      self._domains = [d]
    else:
      raise ValueError("Must specify strand constraints or domain list.")
      
    # Assign length
    self._length = sum([d.length for d in self._domains])
   
   
  ## Basic properties
  @property
  def length(self):
    """ Returns the length of this strand, which equals the sum of the
    lengths of the component domains."""
    return self._length
    
  @property
  def constraints( self ):
    """
    The sequence constraints associated with this strand, computing it as
    necessary from the domains.
    """
    return "".join([d._constraints for d in self.base_domains()])
  @constraints.setter
  def constraints( self, new_constraints ):
    """
    Set the sequence constraints on this Strand. If repeated domains
    are assigned different constraints, the result is the intersection
    of all constraints assigned to the domain.
    """
    assert len(new_constraints) == self._length
    for d in self.base_domains():
      d._constraints = Constraints("N" * d._length)
    self.add_constraints(new_constraints)
  def add_constraints(self, new_constraints):
    """
    Applies the given constraints on top of any existing constraints on the
    sequence of this strand.
    """
    assert len(new_constraints) == self._length
    i = 0
    for d in self.base_domains():
      subconstraints = new_constraints[i : i+d._length]
      d._constraints = d._constraints.intersection(subconstraints)
      i += d._length
  
   
  ## Complementarity and equivalence
  @property
  def is_complement(self):
    """This returns True if the Strand name ends in an asterisk (*).
    False otherwise. ComplementaryStrands negate this."""
    return self.name[-1] == '*'
  @property
  def complement(self):
    """ Returns a ComplementaryStrand object defined by this Strand. """
    return ComplementaryStrand(self)
  def equivalent_to(self, other):
    """ Returns True iff the operand is a strand with the same
    base domain breakdown. """
    return (self._object_type == other._object_type and
            self.base_domains() == other.base_domains())
  def complementary_to(self, other):
    """ Returns True iff the complement of the operand is a strand with
    the same domain breakdown. """
    return (self._object_type == other._object_type and
            self.base_domains() == other.complement.base_domains())
  
  
  ## DNA object hierarchy
  @property
  def domains(self):
    """ Returns the domains originally used to define this strand. 
    If sequence constraints were used to define the strand, returns
    a list with the new domain that holds these constraints."""
    return self._domains
  def base_domains(self):
    """ Returns the unique list of non-composite domains that compose
    this strand."""
    return sum([d.base_domains() for d in self._domains], [])

  ## (In)equality
  def __eq__(self, other):
    """ Two strands are equal if they have the same name. """
    return self._object_type == other._object_type and self.name == other.name
  def __ne__(self, other):
    """ Returns True iff the two strands are not equal."""
    return not self.__eq__(other)
  def __hash__(self):
    """ Returns a hash value on Strands."""
    return sum([ord(c) * pow(2,i)
      for i, c
      in enumerate(self._object_type + self.name)])
  
  ## Output
  def __str__(self):
    """ Human-readable output formatting for this Strand object."""
    info = "[" + ", ".join([str(d) for d in self._domains]) + "]"
    return "Strand {0}: {1} ({2})".format(self.name, info, self._length)
  def __repr__(self):
    return str(self)

class ComplementaryStrand( Strand ):
  """
  Represents a complemented strand. This is always defined in
  terms of an original strand, so that it reflects any change in the
  original. It provides the same interfaces as a strand.
  """

  def __init__( self, complemented_strand ):
    """ Creates a new ComplementaryStrand object defined in
    terms of the given strand."""
    self.id = complemented_strand.id
    
    self._object_type = 'strand'
    
    self._complement = complemented_strand

  ## Basic properties
  @property
  def name( self ):
    """ Returns the name string for this ComplementaryStrand. The name
    is the same as that of its complement, with the presence of the
    trailing asterisk (*) toggled."""
    if self._complement.name.endswith("*") or \
       self._complement.name.endswith("'"):
      return self._complement.name.rstrip("*'")
    else:
      return self._complement.name + "*"
  
  @property
  def length(self):
    """ Returns the length of this ComplementaryStrand, which equals
    the length of its complement."""
    return self._complement.length
    
  @property
  def constraints(self):
    """ Returns the Constraints object for this ComplementaryStrand. This
    is the complementary-reverse of the constraints on its complement."""
    return self._complement.constraints.complement
  @constraints.setter
  def constraints(self, new_constraints):
    """ Sets the Constraints on this ComplementaryStrand by setting those
    of its complement to the complement of the input constraints."""
    self._complement.constraints = new_constraints.complement
  def add_constraints(self, new_constraints):
    """ Applies the given constraints on top of this ComplementaryStrand by
    applying the complement of these constraints onto its complement. """
    self._complement.add_constraints(new_constraints.complement)


  ## Complementarity and equivalence
  @property
  def is_complement(self):
    """ Returns True iff its complement returns False."""
    return not self._complement.is_complement
  @property
  def complement(self):
    """ Returns the complementary strand used to define this
    ComplementaryStrand."""
    return self._complement
  def equivalent_to(self, other):
    """ Returns True iff the RHS is complementary to the LHS's
    complement."""
    return self._complement.complementary_to(other)
  def complementary_to(self, other):
    """ Returns True iff the RHS is equivalent to the LHS's complement."""
    return self._complement.equivalent_to(other)

  ## DNA object hierarchy
  @property
  def domains(self):
    """ Returns the complementary-reverse of the domains that define this
    ComplementaryStrand's complement. """
    return [d.complement for d in reversed(self._complement.domains)]
  def base_domains(self):
    """ Returns the unique list of non-composite domains for this
    ComplementaryStrand. This is complementary-reverse of the list for
    this strand's complement."""
    return [d.complement for d in reversed(self._complement.base_domains())]
    
  ## (In)equality
  def __eq__(self, other):
    """ Returns True iff the complements of the LHS and RHS are equal."""
    return self._complement.__eq__(other.complement)
  def __ne__(self, other):
    """ Returns True iff the LHS and RHS are not equal."""
    return not self.__eq__(other)
  def __hash__(self):
    """ Returns a hash value on ComplementaryStrands. """
    return -self._complement.__hash__()
    
  ## Output
  def __str__(self):
    """ Human-readable output formatting for this ComplementaryStrand."""
    return "ComplementaryStrand {0}: ~[{1}]".format(self.name, str(self._complement))
  def __repr__(self):
    return str(self)
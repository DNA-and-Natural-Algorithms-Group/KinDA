"""
domain.py

Copyright (c) 2010-2014 Caltech. All rights reserved.
Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)
          Joseph Berleant (jberleant@dna.caltech.edu)

This module defines a simple dna `domain` object.

`domain`
    A logical grouping of nucleic acid bases that occur in a specific order. 
"""


from sequence import Sequence


class Domain(object):
  """
  Represents a DNA domain, a sequence of specified length with constraint sequence
  on the allowed nucleotides at each position in the domain.
  """
  id_counter = 0

  def __init__(self, *args, **kargs ):
    """
    Initialization:

    Keyword Arguments:
    name [type=str]                   -- Name of this domain. If not given,
                                         an automatic one is generated.
    sequence [type=str]               -- Sequence constraints on this domain.
                                         Specifiers may be any of
                                         A,T,C,G,R,Y,W,S,M,K,B,V,D,H,N.
                                         Either sequence or subdomains
                                         must be specified.
    subdomains [type=list of Domains] -- List of component Domain objects. Either
                                         sequence or subdomains must be specified.
    """
    # Assign id
    self.id = Domain.id_counter
    Domain.id_counter += 1
    
    # Assign DNA object type
    self._object_type = 'domain'
    
    # Assign name
    if 'name' in kargs: self.name = kargs['name']
    else: self.name = 'domain_{0}'.format(self.id)
    
    # Assign sequence or subdomain list
    if 'sequence' in kargs and 'subdomains' not in kargs:
      self._sequence = Sequence(kargs['sequence'])
      self._length = len(self._sequence)
      self.is_composite = False
    elif 'subdomains' in kargs and 'sequence' not in kargs:
      self._subdomains = kargs['subdomains']
      self._length = sum([d.length for d in self._subdomains])
      self.is_composite = True
    else:
      raise ValueError("Must pass one of 'sequence' or 'subdomains' keyword argument.")
      
      
  ## Basic properties
  @property
  def length(self):
    """ The number of nucleotides in this domain. """
    return self._length
    
  @property
  def sequence(self):
    """Returns the Constraints object or concatenation of the constraint
    string of the subdomains."""
    return Sequence("".join([d._sequence for d in self.base_domains()]))
  @sequence.setter
  def sequence(self, new_seq):
    """Sets the sequence string or sets the sequence of the subdomains.
    If repeated domains are assigned different sequences, the result is
    the intersection of those constraints."""
    assert len(new_seq) == self._length, "KinDA: ERROR: Cannot change domain length after initialization"
    for d in self.base_domains():
      d._sequence = Sequence("N" * d._length)
    self.add_sequence(new_seq)
  def restrict_sequence(self, constraints):
    """Applies the sequence constraints on top of existing sequence constraints on this domain."""
    assert len(constraints) == self._length
    i = 0
    for d in self.base_domains():
      subsequence = constraints[i : i+d._length]
      d._sequence = d._sequence.intersection(subsequence)
      i += d._length

      
  ## Equivalence and complementarity
  @property
  def is_complement(self):
    """This returns True if the Domain name ends in an asterisk (*).
    False otherwise. ComplementaryDomains negate this."""
    return self.name[-1] == '*'
  @property
  def complement(self):
    """
    The complementary domain, specifically the one whose bases in
    the standard 5' to 3' ordering are exactly matching base pairs to
    the original.
    """
    return ComplementaryDomain( self )
  def equivalent_to(self, other):
    """Two domains are equivalent if their base domain lists are equal.
    Note that this is a more generous definition than equality."""
    seqs1 = self.base_domains()
    seqs2 = other.base_domains()
    return seqs1 == seqs2
  def complementary_to(self, other):
    """Two domains are complementary if their base domain lists are the
    complementary-reverse of each other."""
    seqs1 = self.base_domains()
    seqs2_rc = [d.complement for d in reversed(other.base_domains())]
    return seqs1 == seqs2_rc
    

  ## Domain hierarchy
  @property
  def subdomains(self):
    """Returns the subdomains that originally defined this domain. Note
    that these subdomains may be formed from subdomains as well. For
    the maximally split representation of the domain, use base_domains().
    If this domain was originally defined as non-composite, this returns
    a list with only this domain."""
    if self.is_composite:
      return self._subdomains
    else:
      return[self]
  @subdomains.setter
  def subdomains(self, domains):
    """ Sets the subdomains property of this Domain object and sets its
    is_composite flag to True. """
    self._subdomains = domains
    self.is_composite = True
  def base_domains(self):
    """Breaks down the domain into non-composite domains."""
    if self.is_composite:
      return sum([d.base_domains() for d in self._subdomains], [])
    else:
      return [self]
     

  ## (In)equality
  def __eq__(self, other):
    """Returns True iff the two objects are domains and have the same name."""
    return self._object_type == other._object_type and self.name == other.name
  def __ne__(self, other):
    """Returns True iff the two objects are not equal."""
    return not self.__eq__(other)
  def __hash__(self):
    """ Hash function for Domains. """
    return sum([ord(c) * pow(2,i)
      for i, c
      in enumerate(self._object_type + self.name)])
    
  ## Output
  def __str__( self ):
    """Human-readable output formatting for a domain."""
    if self.is_composite:
      info = "[" + ", ".join([str(d.name) for d in self._subdomains]) + "]"
    else:
      info = self._sequence
    return "Domain {0}: {1}".format(self.name, info)
  def __repr__(self):
    return str(self)


class ComplementaryDomain( Domain):
  """
  Represents a complemented domain. Note that this is always
  defined in terms of an original domain and does not have the same
  data members, instead providing an interface to the complementary
  members.

  """
  def __init__(self, complemented_domain ):
    """Create a new ComplementaryDomain in terms of a given domain.
    Note that very little data is stored directly this domain. Most
    data is calculated as requested in terms of the given domain."""
    self._object_type = 'domain'
    self._complement = complemented_domain

  ## Basic properties  
  @property
  def id(self):
    """ The id of a ComplementaryDomain is the id of its complement."""
    return self._complement.id
  @property
  def name(self):
    """ The name of a ComplementaryDomain is the name of its complement with
    the presence of the trailing asterisk (*) toggled. """
    if self._complement.name.endswith("*") or \
       self._complement.name.endswith("'"):
      return self._complement.name.rstrip("*'")
    else:
      return self._complement.name + "*"
  
  @property
  def is_composite(self):
    """ A ComplementaryDomain is composite iff its complement is composite. """
    return self._complement.is_composite
      
  @property
  def length(self):
    """ Returns the length of the ComplementaryDomain, which equals that
    of its complement."""
    return self._complement.length

  @property
  def _sequence(self):
    """ Returns the sequence constraints of this ComplementaryDomain.
    The constraints are the complementary-reverse of the constraints of
    this domain's complement."""
    return self._complement.sequence.complement
  @property
  def sequence(self):
    """ Returns the sequence constraints of this ComplementaryDomain.
    The constraints are the complementary-reverse of the constraints of
    this domain's complement."""
    return self._sequence
  @sequence.setter
  def sequence(self, new_seq):
    """ Sets the sequence constraints of this ComplementaryDomain
    by setting the sequence constraints of its complement to the
    complement of the given constraints."""
    self._complement.sequence = new_seq.complement
  def restrict_sequence(self, constraints):
    """ Applies the given constraints on top of any existing sequence
    constraints on this ComplementaryDomain. """
    self._complement.add_sequence(constraints.complement)
    
      
  ## Equivalence and complementarity
  @property
  def is_complement(self):
    """ Equivalent to not self.complement.is_complement. """
    return not self._complement.is_complement
  @property
  def complement(self):
    """ Returns the domain used to instantiate this ComplementaryDomain. """
    return self._complement
  def equivalent_to(self, other):
    """ Returns true if the RHS is complementary to the complement of the LHS. """
    return self._complement.complementary_to(other)
  def complementary_to(self, other):
    """ Returns true if the RHS is equivalent to the complement of the LHS. """
    return self._complement.equivalent_to(other)
    
    
  ## Domain hierarchy
  @property
  def subdomains(self):
    """ Returns the subdomain list for this ComplementaryDomain. This is
    the complementary-reverse of the subdomain list of its complement."""
    return [d.complement for d in reversed(self._complement.subdomains)]
  @subdomains.setter
  def subdomains(self, domains):
    """ Sets the subdomain list for this ComplementaryDomain and sets its
    is_composite field to True. """
    self._complement.subdomains = [d.complement for d in reversed(domains)]
  def base_domains(self):
    """ Returns the list of non-composite domains that compose this
    ComplementaryDomain. This is the complementary-reverse of the list
    for its complement."""
    return [d.complement for d in reversed(self._complement.base_domains())]

  ## (In)equality
  def __eq__(self, other):
    """ Returns True iff their complements are equal."""
    return self._complement.__eq__(other.complement)
  def __ne__(self, other):
    """ Returns True iff they are not equal."""
    return not self.__eq__(other)
  def __hash__(self):
    """ Hash function for ComplementaryDomains. """
    return -self._complement.__hash__()
    
  ## Output
  def __str__( self ):
    """ Human-readable output formatting for this ComplementaryDomain."""
    if self.is_composite:
      info = "[" + ",".join([d.name for d in self.subdomains]) + "]"
    else:
      info = self.sequence
    return "ComplementaryDomain {0}: {1}".format(self.name, info)
  def __repr__(self):
    return str(self)

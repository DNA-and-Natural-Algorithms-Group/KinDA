"""
domain.py

Copyright (c) 2010 Caltech. All rights reserved.
Coded by: Joseph Schaeffer (schaeffer@dna.caltech.edu)

This module defines a simple dna `domain` object.

`domain`
    A logical grouping of nucleic acid bases that occur in a specific order. 
"""


from strand import Strand

def generate_sequence( n, allowed_bases = ['G','C','T','A'], base_probability = None ):
    """ Generate a sequence of N base pairs.

    Bases are chosen from the allowed_bases [any sequence type], and
    according to the probability distribution base_probability - if
    none is specified, uses uniform distribution."""

    import random
    
    result = ""
    if base_probability is None:
        return result.join([random.choice( allowed_bases ) for i in range(n)])
    else:
        def uniform_seq( r ):
            """ This function returns a lambda to be used on a sequence of tuples of
            (p,item) e.g. via reduce(uniform_seq(.75), seq, (0.0,'none')).

            In this situation, the lambda computes (p',item'), where
            p' is the sum of consecutive probabilities in the sequence
            <= r, and item' is the corresponding item.

            It achieves this by updating the result with the new sum
            of probabilities and last item checked until the the
            summed probability has exceeded the target.

            Note: r should be in [0.0,1.0) as produced by random(),
            any input r >= 1.0 (assuming the sequence given sums to
            1.0) will give the final element as the result which is technically incorrect.
            """
            return lambda x,y: (x[0]+y[0],y[1]) if r>=x[0] else (x[0],x[1])
        return result.join( [reduce( uniform_seq( random.random() ),
                              zip(base_probability,allowed_bases),
                              (0.0,'Invalid Probabilities'))[1] 
                              # note this subscript [1] pulls out the item
                              # selected by the reduce since the result was a tuple.
                      for i in range(n)]
                     )


class Domain(object):
  """
  Represents a sequence domain, for use in defining strands and stop conditions for Multistrand.
  """
  _domain_unique_id = 0

  def __init__(self, *args, **kargs ):
    """
    Initialization:

    Keyword Arguments:
    name [type=str,required]       -- Name of this domain.
    sequence [type=str,default=""] -- Sequence of this domain, e.g. 'AGGCAGATA'
    length [default=0]             -- Length of this domain. A 0-length domain
                                      may be useful in some rare cases. If
                                      sequence is set, length should ALWAYS be
                                      == len(sequence), but this is not strictly
                                      enforced.
    """
    self.id = Domain._domain_unique_id
    Domain._domain_unique_id += 1
    self.sequence = ""
    self.length = 0
    for k,v in kargs.iteritems():
      self.__setattr__( k, v )
    if 'name' not in kargs:
      raise ValueError("Must pass a 'name' keyword argument.")

  def gen_sequence( self, *args, **kargs ):
    """ Uses the same parameters as 'multistrand.utils.generate_sequence', but sets the length to the domain's length."""
    if 'n' in kargs:
      del kargs['n']
    self.sequence = generate_sequence(n = self.length, *args, **kargs )

  @property
  def C(self):
    """
    The complementary domain, specifically the one whose bases in
    the standard 5' to 3' ordering are exactly matching base pairs to
    the original.

    >>> a = Domain( sequence='AGGACCATT')
    >>> a.sequence
    AGGACCATT
    >>> b = a.C
    >>> b.sequence
    AATCCTCCT
    >>> b.name
    >>> b.C == a
    True
    """
    return ComplementaryDomain( self )

  def __add__( self, other ):
    if isinstance(other,Domain):
      return Strand( domains = [self,other] )
    else:
      return NotImplemented

  def __str__( self ):
    return ("\
Domain : {fieldnames[0]:>9}: '{0.name}'\n\
       : {fieldnames[1]:>9}: {0.length}\n".format( self, fieldnames=['Name','Length'] ) +
  (self.sequence and
     '       : {fieldname:>9}: {0}\n'.format( self.sequence, fieldname='Seq' )
     or ''))



class ComplementaryDomain( Domain):
  """
  Represents a complemented domain. Note that this is always
  defined in terms of an original domain and does not have the same
  data members, instead providing an interface to the complementary
  members.

  """
  complement = {'G':'C',
                'C':'G',
                'A':'T',
                'T':'A'}
    
  def __init__(self, complemented_domain ):
    self.id = complemented_domain.id
    
    self._domain = complemented_domain

  @property
  def length(self):
    return self._domain.length

  @property
  def name(self):
    if self._domain.name.endswith("*") or \
       self._domain.name.endswith("'"):
      return self._domain.name.rstrip("*'")
    else:
      return self._domain.name + "*"

  @property
  def sequence( self ):
    if self._domain.sequence == None:
      raise ValueError
    else:
      return "".join([ComplementaryDomain.complement[i] for i in reversed(self._domain.sequence.upper())])

  def gen_sequence( self, *args, **kargs ):
    """ Uses the same parameters as 'multistrand.utils.generate_sequence', but sets the length to the domain's length."""
    self._domain.gen_sequence( *args, **kargs )

  @property
  def C(self):
    """
    The complementary domain, specifically the one whose bases in
    the standard 5' to 3' ordering are exactly matching base pairs to
    the original.

    >>> a = Domain( sequence='AGGACCATT')
    >>> a.sequence
    AGGACCATT
    >>> b = a.C
    >>> b.sequence
    AATCCTCCT
    >>> b.name
    >>> b.C == a
    True
    """
    return self._domain


"""
macrostate.py

Copyright (c) 2014 Caltech. All rights reserved.
Coded by: Joseph Berleant (jberleant@dna.caltech.edu)

This module defines a dna `macrostate`, with structure inspired by
the Multistrand macrostate. Currently, macrostates may only be
formed as a (nested) conjunction or disjunction of five base macrostate types
(see Joseph Schaeffer's PhD thesis on Multistrand for more details):
  EXACT_MACROSTATE        Given a complex of specific structure, this
                          corresponds to all system microstates in which
                          this exact complex exists.
  ORDERED-COMPLEX_MACROSTATE     Given a complex, this corresponds to all system
                          microstates in which there exists a complex with
                          the same strands in the same order.
                          Multistrand calls these DISASSOC macrostates.
  BOUND_MACROSTATE        Given a complex of a single strand, this corresponds
                          to all system microstates in which there exists a
                          complex with this strand and at least one other
                          strand.
  COUNT_MACROSTATE        Given a complex and an integer or percentage, this
                          corresponds to all system microstates in which there
                          exists a complex with a defect less than or equal to
                          this number or percentage of nucleotides in the
                          complex.
  LOOSE_MACROSTATE        Given a complex, an integer or percentage, and a
                          set of nucleotides of interest in the complex, this
                          corresponds to all system microstates in which there
                          is a complex of the same strands and ordering that
                          has a defect over the nucleotides of interest less
                          than or equal to the given number or percentage of
                          nucleotides in the area of interest.

`macrostate`
    Any set of dna microstates (sets of complexes with completely specified structure).
    
"""


class Macrostate(object):
  """ Represents a macrostate as defined by a conjunction of one or
  more macrostates of the above types. """
  
  id_counter = 0
  types = {'exact': 0,
           'bound': 1,
           'ordered-complex': 2,
           'loose': 3,
           'count': 4,
           'conjunction': 5,
           'disjunction': 6}
  
  def __init__(self, *args, **kargs):
    """
    Initialization:
    
    Keyword Arguments:
    name [type=str]                        -- Name of the macrostate. An
                                              automatic name is given if
                                              this is ommitted.
    type [type=str]                        -- Macrostate type. One of
                                              'exact', 'ordered-complex', 'bound',
                                              'count', 'loose', 
                                              'conjunction', or 'disjunction'.
                                              The first 5
                                              are described in Schaeffer's
                                              PhD thesis. The conjunction type
                                              represents a conjunction of
                                              a list of Macrostates. Disjunction is
                                              defined analogously.
    complex [type=str]                     -- Required for all but the
                                              'conjunction' and 'disjunction' macrostates.
    cutoff [type=int OR float]             -- Allowed fractional defect. Required for
                                              'count' and 'loose' macrostates.
    macrostates [type=list of Macrostates] -- Required for the 'conjunction' and 'disjunction'
                                              macrostates.
    """
    # Assign id
    self.id = Macrostate.id_counter
    Macrostate.id_counter += 1
    
    self._object_type = 'macrostate'
    
    # Assign name
    if 'name' in kargs: self.name = kargs['name']
    else: self.name = 'macrostate_{0}'.format(self.id)
    
    # Assign data
    self._macrostate_type = Macrostate.types[kargs['type']]
    if self._macrostate_type == Macrostate.types['exact']:
      self._complex = kargs['complex']
    elif self._macrostate_type == Macrostate.types['ordered-complex']:
      self._complex = kargs['complex']
    elif self._macrostate_type == Macrostate.types['bound']:
      self._complex = kargs['complex']
      assert len(self._complex.strands) == 1, "Bound macrostate requires a single-stranded complex."
    elif self._macrostate_type == Macrostate.types['count']:
      self._complex = kargs['complex']
      self._cutoff = kargs['cutoff']
    elif self._macrostate_type == Macrostate.types['loose']:
      self._complex = kargs['complex']
      self._cutoff = kargs['cutoff']
    elif self._macrostate_type == Macrostate.types['conjunction']:
      self._macrostates = kargs['macrostates']
    elif self._macrostate_type == Macrostate.types['disjunction']:
      self._macrostates = kargs['macrostates']
  
  @property
  def type(self):
    return self._macrostate_type
    
  @property
  def complex(self):
    if self._macrostate_type != Macrostate.types['conjunction'] and self._macrostate_type != Macrostate.types['disjunction']:
      return self._complex
    else:
      raise ValueError("Macrostate lacks a 'complex' field.")
      
  @property
  def cutoff(self):
    if (self._macrostate_type == Macrostate.types['count'] or
        self._macrostate_type == Macrostate.types['loose']):
      return self._cutoff
    else:
      raise ValueError("Macrostate lacks a 'cutoff' field.")
      
  @property
  def macrostates(self):
    if self._macrostate_type == Macrostate.types['conjunction'] or self._macrostate_type == Macrostate.types['disjunction']:
      return self._macrostates
    else:
      raise ValueError("Macrostate lacks a 'macrostates' field.")

  def __str__(self):
    if self._macrostate_type == Macrostate.types['exact']:
      return "Macrostate({}, {})".format('EXACT', str(self._complex))
    elif self._macrostate_type == Macrostate.types['ordered-complex']:
      return "Macrostate({}, {})".format('ORDERED-COMPLEX', str(self._complex))
    elif self._macrostate_type == Macrostate.types['bound']:
      return "Macrostate({}, {})".format('BOUND', str(self._complex))
    elif self._macrostate_type == Macrostate.types['count']:
      return "Macrostate({}, {}, {})".format('COUNT', str(self._complex), self._cutoff)
    elif self._macrostate_type == Macrostate.types['loose']:
      return "Macrostate({}, {}, {})".format('LOOSE', str(self._complex), self._cutoff)
    elif self._macrostate_type == Macrostate.types['conjunction']:
      return "Macrostate({}, {})".format('CONJUNCTION', self._macrostates)
    elif self._macrostate_type == Macrostate.types['disjunction']:
      return "Macrostate({}, {})".format('DISJUNCTION', self._macrostates)
    else:
      return "Macrostate({}, ???)".format(self._macrostate_type)

# kinda.py
# Created by Joseph Berleant, 1/12/2018
#
# Defines the System class, encapsulating statistics calculations for DNA
# strand-displacement system properties.

from typing import Optional

from . import options
from .objects import Complex, RestingSet, Reaction, RestingSetReaction
from .statistics import stats_utils
from .statistics.stats import RestingSetStats, RestingSetRxnStats


class System:
  """
  Stores and manages Stats objects for each system component for easier
  retrieval. Data can also be stored in a file and retrieved for later analysis.
  A System object is instantiated with an EnumerateJob object, from which
  detailed and condensed reactions as well as resting sets and complexes are
  taken.
  """
  def __init__(self,
               complexes, 
               restingsets = [], 
               detailed_reactions = [], 
               condensed_reactions = [], 
               enumeration = True, 
               kinda_params = {}, 
               peppercorn_params = {}, 
               multistrand_params = {}, 
               nupack_params = {}):
    """
    Constructs a System object with the given complexes, restingsets, reactions,
    and condensed reactions. If enumeration is True (default), the Peppercorn
    enumerator is used to enumerate a detailed reaction network. System performs
    reaction condensation to produce the condensed reaction network. The given
    restingsets and reactions are added to the enumerated network. Detailed
    reactions are added to the network prior to condensation. If enumeration is
    False, then no enumeration is performed and only the given resting sets and
    reactions are used for analysis.
    """
    # Update parameters for KinDA, Multistrand, and Nupack
    self._kinda_params = dict(options.kinda_params, **kinda_params)
    self._multistrand_params = dict(options.multistrand_params, **multistrand_params)
    self._nupack_params = dict(options.nupack_params, **nupack_params)

    r_model, r_method, uni_sc, bi_sc = [
      self._multistrand_params.get(k, None) for k in
      ["rate_model", "rate_method", "unimolecular_scaling", "bimolecular_scaling"]]
    assert (
      (r_model is None and isinstance(r_method, str)
       and all(isinstance(k, float) for k in [uni_sc, bi_sc])) or
      (isinstance(r_model, str)
       and all(k is None for k in [r_method, uni_sc, bi_sc]))), (
        "In `multistrand_params`, please specify either one of:\n"
        "  - rate_model: str\n"
        "  - rate_method: str, "
        "unimolecular_scaling: float, bimolecular_scaling: float.")

    # Store initial DSD system objects
    assert all(isinstance(rxn, RestingSetReaction) for rxn in condensed_reactions)
    self._condensed_reactions = set(condensed_reactions)
    assert all(isinstance(rxn, Reaction) for rxn in detailed_reactions)
    self._detailed_reactions = set(detailed_reactions)

    assert all(isinstance(rs, RestingSet) for rs in restingsets)
    self._restingsets = set(restingsets
        + [rs for rxn in self._condensed_reactions
           for rs in rxn.reactants+rxn.products])
    assert all(isinstance(c, Complex) for c in complexes)
    self._complexes = set(complexes
        + [c for rs in self._restingsets for c in rs.complexes]
        + [c for rxn in self._detailed_reactions
           for c in rxn.reactants+rxn.products])

    if enumeration:
      self._peppercorn_params = dict(options.peppercorn_params, **peppercorn_params)
      self.enumerate()
    else:
      self._peppercorn_params = {}
      self._enum_job = None

    # Make stats objects separately, not during initialization. You may want to 
    # specify or query parameters before ...
    self._rs_to_stats = None
    self._rxn_to_stats = None
    self._spurious_restingsets = None
    self._spurious_condensed_reactions = None

    # ok, I changed my mind, let's make stats objects right away ...  but I can
    # still see how one wants to initialize the object, then twiggle some
    # session parameters and *then* make the stats_objects.
    self.make_stats_objects()

  def enumerate(self):
    from .enumeration.enumeratejob import EnumerateJob

    # Create enumeration object
    enum_job = EnumerateJob(
        complexes = self.complexes,
        reactions = self.detailed_reactions,
        peppercorn_params = self._peppercorn_params
    )

    # Incorporate enumerated data
    self._enum_job = enum_job
    self._complexes = set(enum_job.get_complexes())
    self._restingsets = set(enum_job.get_restingsets())
    self._detailed_reactions = set(enum_job.get_reactions())
    self._condensed_reactions = set(enum_job.get_restingset_reactions())

  def make_stats_objects(self):
    """
    Create stats objects for reactions and resting sets.
    `make_stats()` will also make stats objects for potential spurious reactions
    and resting sets.
    """
    # Filter out unimolecular reactions, if these are disabled.
    if not self._kinda_params['enable_unimolecular_reactions']:
      self._condensed_reactions = set([
        rxn for rxn in self._condensed_reactions if len(rxn.reactants) == 2])

    # Create stats objects for reactions and resting sets make_stats() will also
    # make stats objects for potential spurious reactions and resting sets
    self._rs_to_stats, self._rxn_to_stats = stats_utils.make_stats(
        list(self._complexes),
        list(self._restingsets),
        list(self._detailed_reactions),
        list(self._condensed_reactions),
        kinda_params = self._kinda_params,
        multistrand_params = self._multistrand_params,
        nupack_params = self._nupack_params
    )
    # Pull out spurious reactions/resting sets
    self._spurious_restingsets = set(self._rs_to_stats.keys()) - self._restingsets
    self._spurious_condensed_reactions = (
      set(self._rxn_to_stats.keys()) - self._condensed_reactions)

    # Set default max concentration for each resting set
    for rs_stats in self._rs_to_stats.values():
      rs_stats.c_max = self._kinda_params['max_concentration']

  ## Basic get functions for system objects
  @property
  def initialization_params(self):
    """
    Returns a dict of the parameters used to initialize this system, including
    those set by the defaults in options.py.
    """
    return {
      'kinda_params': self._kinda_params.copy(),
      'multistrand_params': self._multistrand_params.copy(),
      'nupack_params': self._nupack_params.copy(),
      'peppercorn_params': self._peppercorn_params.copy()
    }

  @property
  def kinda_params(self):
    """
    Returns a dict of the KinDA parameters used when initializing the system.
    Equivalent to initialization_params['kinda_params'].
    """
    return self._kinda_params.copy()

  @property
  def multistrand_params(self):
    """
    Returns a dict of the Multistrand parameters used when initializing the
    system. Equivalent to initialization_params['multistrand_params'].
    """
    return self._multistrand_params.copy()

  @property
  def nupack_params(self):
    """
    Returns a dict of the NUPACK parameters used when initializing the system.
    Equivalent to initialization_params['nupack_params'].
    """
    return self._nupack_params.copy()

  @property
  def peppercorn_params(self):
    """
    Returns a dict of the Peppercorn parameters used when initializing the
    system. Equivalent to the initialization_params['peppercorn_params'].
    """
    return self._peppercorn_params.copy()

  @property
  def complexes(self):
    """
    Returns a list of all complexes (given and enumerated) predicted for the
    system.
    """
    return list(self._complexes)

  @property
  def restingsets(self):
    """ Returns a list of all resting sets (given, enumerated, and spurious)
    predicted for the system. """
    return list(self._restingsets | self._spurious_restingsets)

  @property
  def detailed_reactions(self):
    """
    Returns a list of detailed reactions (given and enumerated) predicted for
    the system.
    """
    return list(self._detailed_reactions)

  @property
  def condensed_reactions(self):
    """
    Returns a list of all resting-set (condensed) reactions
    (given, enumerated, and spurious) for the system.
    """
    return list(self._condensed_reactions | self._spurious_condensed_reactions)

  ## Convenience filters for specific objects
  def get_reactions(self, reactants = [], products = [], arity = 2,
                    unproductive = None, spurious = None):
    """
    Returns a list of all reactions including the given reactants and the given
    products. If specified, spurious = True will return only spurious reactions
    (those not enumerated by Peppercorn) and spurious = False will return only
    enumerated reactions. Otherwise, no distinction will be made.
    """
    if spurious == True:
      rxns = list(self._spurious_condensed_reactions)
    elif spurious == False:
      rxns = list(self._condensed_reactions)
    else:
      rxns = list(self._spurious_condensed_reactions | self._condensed_reactions)

    if unproductive == True:
      rxns = [x for x in rxns if x.has_reactants(x.products) and x.has_products(x.reactants)]
    elif unproductive == False:
      rxns = [x for x in rxns if not(x.has_reactants(x.products) and x.has_products(x.reactants))]

    if arity is not None:
      rxns = [x for x in rxns if len(x.reactants)==arity]

    return sorted([x for x in rxns
                   if x.has_reactants(reactants) and x.has_products(products)])

  def get_reaction(self, **kwargs):
    """ Returns a single reaction matching the criteria given. """
    rxns = self.get_reactions(**kwargs)
    if len(rxns) == 0:
      print("KinDA: ERROR: SystemStats.get_reaction() failed to find "
            "a reaction with the given criteria.")
      return None
    elif len(rxns) > 1:
      print("KinDA: WARNING: SystemStats.get_reactino() found "
            "multiple reactions with the given criteria.")
    return rxns[0]

  def get_restingsets(self, complex = None, strands = [], name = None,
                      complex_name = None, spurious = False):
    """ Returns a list of resting sets satisfying the filter arguments.

    Args:
        complex (optional): Returns a list of resting sets containing the given complex.
        strands (optional): Returns a list of resting sets containing the given strands.
        name (optional): Returns a list of resting sets with the given name.
        complex_name (optional): Returns a list of resting sets containing a complex with this name.
        suprious (optional): True: Returns only spurious resting sets.
                             False: Returns only non-spurious resting sets.
                             None: Returns both spurious and non-spurious resting sets.
    """
    if spurious == True:
      rs = list(self._spurious_restingsets)
    elif spurious == False:
      rs = list(self._restingsets)
    else:
      rs = list(self._restingsets | self._spurious_restingsets)

    if complex is not None:
      rs = [x for x in rs if complex in x]
    if strands != []:
      rs = [x for x in rs if all([s in x.strands for s in strands])]
    if name is not None:
      rs = [x for x in rs if x.name == name]
    if complex_name is not None:
      rs = [x for x in rs if complex_name in [c.name for c in x.complexes]]

    return sorted(rs)

  def get_restingset(self, complex = None, strands = [], name = None,
                     complex_name = None, spurious = False):
    rs_list = self.get_restingsets(
      complex = complex, strands = strands, name = name,
      complex_name = complex_name, spurious = spurious)
    if len(rs_list) == 0:
      print("KinDA: ERROR: SystemStats.get_restingset() failed to find "
            "a resting set with the given criteria")
      return None
    elif len(rs_list) > 1:
      print("KinDA: WARNING: SystemStats.get_restingset() found "
            "multiple resting sets with the given criteria")
    return rs_list[0]

  def get_complexes(self, name = None):
    complexes = list(self._complexes)
    if name is not None:
      complexes = [x for x in complexes if x.name == name]
    return sorted(complexes)

  def get_complex(self, name = None):
    complexes = self.get_complexes(name = name)
    if len(complexes) == 0:
      print("KinDA: ERROR: SystemStats.get_complexes() failed to find "
            "a complex with the given criteria.")
      return None
    elif len(complexes) > 1:
      print("KinDA: WARNING: SystemStats.get_complexes() found "
            "multiple complexes with the given criteria.")
    return complexes[0]

  def get_stats(self, obj: RestingSet | RestingSetReaction
                ) -> Optional[RestingSetStats | RestingSetRxnStats]:
    """
    Returns the stats object corresponding to the given system object.
    """
    assert self._rs_to_stats is not None
    assert self._rxn_to_stats is not None
    assert isinstance(obj, (RestingSet, RestingSetReaction))

    if isinstance(obj, RestingSet) and obj in self._rs_to_stats:
      return self._rs_to_stats[obj]
    elif isinstance(obj, RestingSetReaction) and obj in self._rxn_to_stats:
      return self._rxn_to_stats[obj]
    else:
      print(f"Statistics for object not found: {obj}")
      return None

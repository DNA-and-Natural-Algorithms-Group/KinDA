### TODO: Add conversion from Multistrand objects to DNAObjects objects

import itertools as it
from . import utils
##################
## Use to_Multistrand() to convert a system of domains, strands, complexes, resting sets, and
## macrostates to equivalent Multistrand objects.
##################

def to_Multistrand_domains(domains):
  import multistrand.objects as MS
  
  ## Pull out just the base domains, which are the ones included in the dict
  base_domains = set(sum([d.base_domains() for d in domains], []))
  ## Remove redundant complementary pairs of complexes (we only need one or the other)
  base_domains -= set([d.complement for d in base_domains if not d.is_complement])

  ## Create dict of Multistrand Domain objects
  ms_domains = {}
  for d in base_domains:
    seq = utils.random_sequence(d.sequence)
    ms_domains[d] = MS.Domain(name = d.name, sequence = seq, length = d.length)
    ms_domains[d.complement] = ms_domains[d].C
  return ms_domains
  
def to_Multistrand_strands(strands, ms_domains):
  import multistrand.objects as MS
  
  ms_strands = {}
  for s in strands:
    if not s.is_complement:
      strand_domains = [ms_domains[d] for d in s.base_domains()]
      ms_strands[s] = MS.Strand(name = s.name, domains = strand_domains)
      ms_strands[s.complement] = ms_strands[s].C
  return ms_strands    
  
def to_Multistrand_complexes(complexes, ms_strands):
  import multistrand.objects as MS
  
  ms_complexes = {}
  for c in complexes:
    complex_strands = [ms_strands[s] for s in c.strands]
    struct = c.structure.to_dotparen()
    ms_complexes[c] = MS.Complex(name = c.name,
                                 strands = complex_strands,
                                 structure = struct)
  return ms_complexes
  
def to_Multistrand_restingstates(resting_sets, ms_complexes):
  """ Because support for Multistrand RestingStates has been dropped,
  we now convert to Multistrand Complex objects with the structure
  of one of the component complexes in the KinDA RestingSet object. """
  import multistrand.objects as MS
  
  ms_restingstates = {}
  for rs in resting_sets:
    ms_restingstates[rs] = ms_complexes[next(iter(rs.complexes))]
  return ms_restingstates
  
def to_Multistrand_macrostates(macrostates, ms_complexes):
  import multistrand.objects as MS
  from .macrostate import Macrostate
  from .utils import num_wildcards, macrostate_to_dnf
  
  EXACT = 0
  BOUND = 1
  DISASSOC = 2
  LOOSE = 3
  COUNT = 4
  ms_macrostates = {}
  for m in macrostates:
    m_dnf = macrostate_to_dnf(m)

    if m_dnf.type == Macrostate.types['exact']:
      c = ms_complexes[m_dnf.complex]
      ms_macrostates[m] = [MS.Macrostate(m.name, [(c, EXACT, 0)])]
    elif m_dnf.type == Macrostate.types['ordered-complex']:
      c = ms_complexes[m_dnf.complex]
      ms_macrostates[m] = [MS.Macrostate(m.name, [(c, DISASSOC, 0)])]
    elif m_dnf.type == Macrostate.types['bound']:
      c = ms_complexes[m_dnf.complex]
      ms_macrostates[m] = [MS.Macrostate(m.name, [(c, BOUND, 0)])]
    elif m_dnf.type == Macrostate.types['count']:
      c = ms_complexes[m_dnf.complex]
      if isinstance(m_dnf.cutoff, float):
        cutoff = int(m_dnf.cutoff * (len(c.structure) - c.structure.count('*')))
      else:
        cutoff = m_dnf.cutoff
      ms_macrostates[m] = [MS.Macrostate(m.name, [(c, COUNT, cutoff)])]
    elif m_dnf.type == Macrostate.types['loose']:
      c = ms_complexes[m_dnf.complex]
      if isinstance(m_dnf.cutoff, float):
        cutoff = int(m_dnf.cutoff * (len(c.structure) - c.structure.count('*')))
      else:
        cutoff = m_dnf.cutoff
      ms_macrostates[m] = [MS.Macrostate(m.name, [(c, LOOSE, cutoff)])]
    elif m_dnf.type == Macrostate.types['conjunction']:
      for ms in m_dnf.macrostates:
        assert(ms.type != Macrostate.types['conjunction'] and ms.type != Macrostate.types['disjunction'])
        if ms not in ms_macrostates:
          ms_macrostates[ms] = list(to_Multistrand_macrostates([ms], ms_complexes).values())[0]
      states = sum([ms_macrostates[s][0].complex_items for s in m_dnf.macrostates], [])
      ms_macrostates[m] = [MS.Macrostate(m.name, states)]
    elif m_dnf.type == Macrostate.types['disjunction']:
      ms_macrostates[m] = list(it.chain(*list(to_Multistrand_macrostates(m_dnf.macrostates, ms_complexes).values())))
      for ms in ms_macrostates[m]:  ms.tag = m.name
  return ms_macrostates
  
def to_Multistrand(*args, **kargs):
  """Converts DNAObjects objects to Peppercorn objects.
  The following objects will be recognized in the key list and converted:
    - domains
    - strands
    - complexes
    - resting_sets
    - macrostates
  Other supplied objects are ignored.
  The output is in the form of a dictionary mapping the keys
  'domains', 'strands', 'complexes', 'resting_sets', and macrostates
  to a list of tuples of the form (object, converted_object). This may be
  iterated through directly or converted into a dict for object lookup.
  """
  import multistrand.objects as MS
  from .macrostate import Macrostate
  from .utils import get_dependent_complexes
  
  ## Extract all objects to be converted
  macrostates = set(kargs.get('macrostates', []))
  resting_sets = set(kargs.get('resting_sets', []))
  complexes = set(sum([get_dependent_complexes(m) for m in macrostates], [])
                  + sum([list(rs.complexes) for rs in resting_sets], [])
                  + kargs.get('complexes', []))
  strands = set(sum([c.strands for c in complexes], [])
                + kargs.get('strands', []))
  domains = set(sum([s.base_domains() for s in strands], [])
                + sum([d.base_domains() for d in kargs.get('domains', [])], []))

  ## Create dict of Multistrand Domain objects
  ms_domains = to_Multistrand_domains(domains)
    
  ## Create dict of Multistrand Strand objects
  ms_strands = to_Multistrand_strands(strands, ms_domains)
      
  ## Create dict of Multistrand Complex objects
  ms_complexes = to_Multistrand_complexes(complexes, ms_strands)
  
  ## Create dict of Multistrand RestingState objects
  ms_restingstates = to_Multistrand_restingstates(resting_sets, ms_complexes)
  
  ## Create dict of Multistrand Macrostate objects
  ms_macrostates = to_Multistrand_macrostates(macrostates, ms_complexes)
    
  results = {'domains': list(ms_domains.items()),
             'strands': list(ms_strands.items()),
             'complexes': list(ms_complexes.items()),
             'restingstates': list(ms_restingstates.items()),
             'macrostates': list(ms_macrostates.items())}
  return results
    
  

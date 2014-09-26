from imports import multistrandhome

import utils
##################
## Use to_Multistrand() to convert a system of domains, strands, complexes, resting sets, and
## macrostates to equivalent Multistrand objects.
##################

def to_Multistrand_domains(domains):
  import multistrand.objects as MS
  
  ## Create dict of Multistrand Domain objects
  ms_domains = {}
  for d in domains:
    if not d.is_composite and not d.is_complement:
      seq = utils.random_sequence(d.constraints)
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
  import multistrand.objects as MS
  
  ms_restingstates = {}
  for rs in resting_sets:
    cpxs = [ms_complexes[c] for c in rs.complexes]
    ms_restingstates[rs] = MS.RestingState(name = rs.name, complex_set = cpxs)
  return ms_restingstates
  
def to_Multistrand_macrostates(macrostates, ms_complexes):
  import multistrand.objects as MS
  from macrostate import Macrostate
  from utils import num_wildcards
  
  EXACT = 0
  BOUND = 1
  DISASSOC = 2
  LOOSE = 3
  COUNT = 4
  ms_macrostates = {}
  for m in macrostates:
    if m.type == Macrostate.types['exact']:
      c = ms_complexes[m.complex]
      ms_macrostates[m] = MS.Macrostate(m.name, [(c, EXACT, 0)])
    elif m.type == Macrostate.types['disassoc']:
      c = ms_complexes[m.complex]
      ms_macrostates[m] = MS.Macrostate(m.name, [(c, DISASSOC, 0)])
    elif m.type == Macrostate.types['bound']:
      c = ms_complexes[m.complex]
      ms_macrostates[m] = MS.Macrostate(m.name, [(c, BOUND, 0)])
    elif m.type == Macrostate.types['count']:
      c = ms_complexes[m.complex]
      if isinstance(m.cutoff, float):
        cutoff = int(m.cutoff * (len(c.structure) - c.structure.count('*')))
      else:
        cutoff = m.cutoff
      ms_macrostates[m] = MS.Macrostate(m.name, [(c, COUNT, cutoff)])
    elif m.type == Macrostate.types['loose']:
      c = ms_complexes[m.complex]
      if isinstance(m.cutoff, float):
        cutoff = int(m.cutoff * (len(c.structure) - c.structure.count('*')))
      else:
        cutoff = m.cutoff
      ms_macrostates[m] = MS.Macrostate(m.name, [(c, LOOSE, cutoff)])
    elif m.type == Macrostate.types['conjunction']:
      for ms in m.macrostates:
        if ms not in ms_macrostates:
          ms_macrostates[ms] = to_Multistrand_macrostates([ms], ms_complexes).values()[0]
      states = sum([ms_macrostates[s].complex_items for s in m.macrostates], [])
      ms_macrostates[m] = MS.Macrostate(m.name, states)
  return ms_macrostates
  
def to_Multistrand(domains, strands = [], complexes = [], resting_sets = [], macrostates = []):
  import multistrand.objects as MS
  from macrostate import Macrostate
  
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
    
  return (ms_domains, ms_strands, ms_complexes, ms_restingstates, ms_macrostates)
    
  
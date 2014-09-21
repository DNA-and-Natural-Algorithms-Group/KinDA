from macrostate import Macrostate

## Sequence utilities
# Global DNA nucelotide groups
base_group =      {"A": "A",   "T": "T",   "C": "C",   "G": "G",
                   "R": "AG",  "Y": "CT",  "W": "AT",  "S": "CG",  "M": "AC",  "K": "GT", 
                   "B": "CGT", "V": "ACG", "D": "AGT", "H": "ACT",
                   "N": "ACGT", "x": ""
                  }
base_complement = {"A": "T",   "T": "A",   "C": "G",   "G": "C",
                   "R": "Y",   "Y": "R",   "W": "W",   "S": "S",   "M": "K",   "K": "M",
                   "B": "V",   "V": "B",   "D": "H",   "H": "D",
                   "N": "N",   "x": "x"
                  }
                  
def random_sequence(constraint, base_probs = None):
  import random
  def choose_base(bases, base_probs):
    prob_sum = sum([base_probs[b] for b in bases])
    assert prob_sum > 0, "Cannot obtain sequence with constraint " + constraint
    
    n = random.uniform(0, prob_sum)
    for b in bases:
      n -= base_probs[b]
      if n <= 0: return b
    
  if base_probs == None:
    base_probs = {'A': 0.2, 'T': 0.2, 'C': 0.3, 'G': 0.3}
    
  return "".join(map(lambda b: choose_base(base_group[b], base_probs), constraint))
    
## Functions for Complex objects
def defect(complex, structure):
  """ Calculates the defect of the complex's structure against the given
  target structure. The defect is the number of nucleotides bound
  differently from the target. """
  defect = 0
  for strand_num, strand in enumerate(complex.strands):
    for i in range(strand.length):
      if (complex.bound_to(strand_num, i) != structure.bound_to(strand_num, i)
          and structure.bound_to(strand_num, i) != '?'):
        defect += 1
  return defect
  
def domain_defect(complex, strand_num, domain_num, structure):
  """ Calculates the defect of the binding of the given domain in
  the complex against the domain's expected binding given in the
  specified target structure. """
  strandlist = structure.structure
  strands = structure.strands
  
  for strand_n, strand_struct in enumerate(strandlist):
    if strand_n != strand_num:
      strandlist[strand_n] = ['?'] * len(strand_struct)
  domain_start = sum([d.length
                        for d
                        in strands[strand_num].domains[:domain_num]])
  domain_end = domain_start + strands.domains[domain_num].length
  strandlist[strand_num][:domain_start] = ['?'] * domain_start
  strandlist[strand_num][domain_end:] = ['?'] * (strands[strand_num].length - domain_end)
  
  new_struct = Structure(structure = strandlist, strands = structure.strands)
  return defect(complex, new_struct)
  
def max_domain_defect(complex, structure):
  """ Returns the maximum domain defect among all domains in the given
  complex when compared against the given structure. The defects are
  calculated as decimal fractions of the length of each domain, to
  make the defects comparable across different domains. """
  domains = complex.base_domains()
  max_defect = 0
  for strand_num, strand_domains in enumerate(domains):
    for domain_num, domain in enumerate(strand_domains):
      if domain.length > 0:
        defect = domain_defect(complex, strand_num, domain_num, structure) / float(domain.length)
        max_defect = max(defect, max_defect)
  return max_defect
  
  
## Functions for Macrostates
def get_dependent_complexes(macrostate):
  """ Returns a list of the Complex objects upon which this
  Macrostate depends. """
  if macrostate.type == Macrostate.types['conjunction']:
    return list(set(sum([get_dependent_complexes(m) for m in macrostate.macrostates], [])))
  else:
    return [macrostate.complex]
  
  
## Multistrand utilities
def to_Multistrand(domains, strands = [], complexes = [], resting_sets = [], macrostates = []):
  import multistrand.objects as MS
  
  ## Create dict of Multistrand Domain objects
  ms_domains = {}
  for d in domains:
    if not d.is_composite and not d.is_complement:
      seq = random_sequence(d.constraints)
      ms_domains[d] = MS.Domain(name = d.name, sequence = seq, length = d.length)
      ms_domains[d.complement] = ms_domains[d].C
    
  ## Create dict of Multistrand Strand objects
  ms_strands = {}
  for s in strands:
    if not s.is_complement:
      strand_domains = [ms_domains[d] for d in s.base_domains()]
      ms_strands[s] = MS.Strand(name = s.name, domains = strand_domains)
      ms_strands[s.complement] = ms_strands[s].C
      
  ## Create dict of Multistrand Complex objects
  ms_complexes = {}
  for c in complexes:
    complex_strands = [ms_strands[s] for s in c.strands]
    struct = c.structure.to_dotparen()
    ms_complexes[c] = MS.Complex(name = c.name,
                                 strands = complex_strands,
                                 structure = struct)
  
  ## Create dict of Multistrand RestingState objects
  ms_restingstates = {}
  for rs in resting_sets:
    cpxs = [ms_complexes[c] for c in rs.complexes]
    ms_restingstates[rs] = MS.RestingState(name = rs.name, complex_set = cpxs)
  
  ## Create dict of Multistrand Macrostate objects
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
      ms_macrostates[m] = MS.Macrostate(m.name, [(c, COUNT, m.cutoff)])
    elif m.type == Macrostate.types['loose']:
      c = ms_complexes[m.complex]
      ms_macrostates[m] = MS.Macrostate(m.name, [(c, LOOSE, m.cutoff)])
    elif m.type == Macrostate.types['conjunction']:
      states = sum([ms_macrostates[s].complex_items for s in m.macrostates], [])
      ms_macrostates[m] = MS.Macrostate(m.name, states)
    
  return (ms_domains, ms_strands, ms_complexes, ms_restingstates, ms_macrostates)
    
  
import itertools as it

## Sequence utilities
                  
def random_sequence(sequence, base_probs = None):
  import random
  from .sequence import base_group
  def choose_base(bases, base_probs):
    prob_sum = sum([base_probs[b] for b in bases])
    assert prob_sum > 0, "Cannot obtain sequence with constraint " + sequence
    
    n = random.uniform(0, prob_sum)
    for b in bases:
      n -= base_probs[b]
      if n <= 0: return b
    
  if base_probs == None:
    base_probs = {'A': 0.2, 'T': 0.2, 'C': 0.3, 'G': 0.3}
    
  return "".join([choose_base(base_group[b], base_probs) for b in sequence])
  
## Functions for Domain objects
def split_domain(domain, pos):  
  """ Splits a domain into two subdomains. """
  from .domain import Domain
  
  section1 = Domain(name = domain.name + "[:{0}]".format(pos),
                    sequence = domain.sequence[:pos])
  section2 = Domain(name = domain.name + "[{0}:]".format(pos + 1),
                    sequence = domain.sequence[pos:])
      
  domain.subdomains = [section1, section2]

def equate_domains(domains):
  """This is pretty slow and inefficient."""
  from .domain import Domain
  
  def update_subdomains_list(subdomains):
    new_subdomains = []
    for i, d in enumerate(subdomains):
      if d.is_composite:
        new_subdomains.extend(d.subdomains)
      else:
        new_subdomains.append(d)
    return new_subdomains
  def equate_coincident_domains(domain1, domain2):
    assert domain1.length == domain2.length
    
    # Make a new Domain with sequence from both domains
    new_domain = Domain(name = "~(" + min(domain1.name, domain2.name) + ")",
                        sequence = domain1.sequence)
    new_domain.restrict_sequence(domain2.sequence)
    
    # Redirect domain1 and domain2 to point to new domain
    domain1.subdomains = [new_domain]
    domain2.subdomains = [new_domain]
  
  domains = list(domains)
  index = 1
  while index < len(domains):
    domains1 = domains[index-1].base_domains()
    domains2 = domains[index].base_domains()
    
    # Iterate through each the sub-domain lists for each Domain object,
    # splitting up larger domain objects as needed
    while len(domains1) != 0:
      d1 = domains1[0]
      d2 = domains2[0]
      if d1.length > d2.length:
        # The first subdomain of the first Domain is longer, so split it
        split_domain(d1, d2.length)
        equate_coincident_domains(d1.subdomains[0], d2)
      elif d2.length > d1.length:
        # The first subdomain of the second Domain is longer, so split it
        split_domain(d2, d1.length)
        equate_coincident_domains(d2.subdomains[0], d1)
      elif d1.complementary_to(d2):
        assert(d1.length % 2 == 0), "Unable to equate odd-length %s with its complement." % d1
        split_domain(d1, d1.length / 2) # Split sequence in half, make the halves complementary
        equate_coincident_domains(d1.subdomains[0], d1.subdomains[1].complement)
      elif d1 != d2: # lengths are equal but they are not the same domain
        # The domains coincide, so make them redirect to a common domain
        equate_coincident_domains(d1, d2)
        
      domains1 = update_subdomains_list(domains1)[1:]
      domains2 = update_subdomains_list(domains2)[1:]
        
    index+=1
  
    
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
  specified target structure. If the domain is not completely bound
  or unbound, the results are not very meaningful.
  The defect is calculated as a fraction of the domain's length."""
  from .structure import Structure
  
  strandlist = [lst[:] for lst in structure.to_strandlist()]
  strands = structure.strands

  domain_start = sum([d.length
                        for d
                        in strands[strand_num].domains[:domain_num]])
  domain_end = domain_start + strands[strand_num].domains[domain_num].length
  n = 0
  for strand_n, strand_struct in enumerate(strandlist):
    for i, bound in enumerate(strand_struct):
      is_domain = (strand_n == strand_num and domain_start <= i < domain_end)
      is_bound_to_domain = (bound != None and bound != '?'
                           and bound[0] == strand_num
                           and domain_start <= bound[1] < domain_end)
      if not is_domain and not is_bound_to_domain:
        strandlist[strand_n][i] = '?'
      else:
        n += 1
        
  new_struct = Structure(structure = strandlist, strands = structure.strands)
  return defect(complex, new_struct) / float(n)
  
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
        defect = domain_defect(complex, strand_num, domain_num, structure)
        max_defect = max(defect, max_defect)
  return max_defect
  
  
## Functions for Macrostates
def get_dependent_complexes(macrostate):
  """ Returns a list of the Complex objects upon which this
  Macrostate depends. """
  from .macrostate import Macrostate
  
  if macrostate.type == Macrostate.types['conjunction'] or macrostate.type == Macrostate.types['disjunction']:
    return list(set(sum([get_dependent_complexes(m) for m in macrostate.macrostates], [])))
  else:
    return [macrostate.complex]
def macrostate_to_dnf(macrostate, simplify = True):
  """ Returns a macrostate in disjunctive normal form (i.e. an OR of ANDs).
  Note that this may lead to exponential explosion in the number of terms.
  However it is necessary when creating Multistrand Macrostates, which can
  only be represented in this way. Also, we don't try to simplify much 
  so the expressions may be inefficient/redundant. Adding simplifications of the
  logical expression using (e.g.) De Morgan's laws is a future optimization. """
  from .macrostate import Macrostate

  if macrostate.type != Macrostate.types['conjunction'] and macrostate.type != Macrostate.types['disjunction']:
    dnf_macrostates = [Macrostate(type='conjunction', macrostates=[macrostate])]
  elif macrostate.type == Macrostate.types['conjunction']:
    clauses = [macrostate_to_dnf(m, simplify=False) for m in macrostate.macrostates]
    dnf_macrostates = clauses[0].macrostates
    for clause in clauses[1:]:
      # multiply two dnf clauses
      dnf_macrostates = [Macrostate(type='conjunction', macrostates=m1.macrostates+m2.macrostates) for m1,m2 in it.product(dnf_macrostates, clause.macrostates)]
  elif macrostate.type == Macrostate.types['disjunction']:
    clauses = [macrostate_to_dnf(m, simplify=False) for m in macrostate.macrostates]
    dnf_macrostates = []
    for clause in clauses:
      # add two dnf clauses
      dnf_macrostates += clause.macrostates

  # The most basic simplification. We just subsitute AND/OR expressions with only one operand
  # with just that operand.
  if simplify:
    for i,m in enumerate(dnf_macrostates):
      if len(m.macrostates) == 1:  dnf_macrostates[i]=m.macrostates[0]
  if simplify and len(dnf_macrostates)==1:
    dnf = dnf_macrostates[0]
  else:
    dnf = Macrostate(type='disjunction', macrostates=dnf_macrostates)

  return dnf
def print_macrostate_tree(m, prefix=''):
  from .macrostate import Macrostate

  print("{}|".format(prefix))
  print("{}|".format(prefix))

  if m.type != Macrostate.types['conjunction'] and m.type != Macrostate.types['disjunction']:
    print("{}{}".format(prefix, str(m)))
    print(prefix)
    return
  
  char = '*' if m.type==Macrostate.types['conjunction'] else '+'
  for i,v in enumerate(m.macrostates):
    print("{}{}-- ".format(prefix, char))
    if i < len(m.macrostates)-1:
      print_macrostate_tree(v, prefix+"|  ")
    else:
      print_macrostate_tree(v, prefix+"   ")


    
## Complex -> Macrostate functions
def exact_complex_macrostate(complex):
  from .macrostate import Macrostate

  return Macrostate(name = "macrostate_" + complex.name,
                    type = 'exact',
                    complex = complex)
def count_by_complex_macrostate(complex, cutoff):
  """ Creates a macrostate corresponding to secondary structures that
  match the binding of the given complex, within the given cutoff.
  cutoff is a fractional defect over the entire complex. """
  from .macrostate import Macrostate

  return Macrostate(name = "macrostate_{}_{}".format(complex.name, cutoff),
                    type = "count",
                    complex = complex,
                    cutoff = cutoff)
def loose_domain_macrostate(complex, strand_num, domain_num, cutoff):
  """ Creates a loose macrostate with the given domain as the region
  of interest. The cutoff is a fractional defect within this domain. The
  results are not meaningful if the domain is not completely bound or 
  completely unbound."""
  from .complex import Complex
  from .structure import Structure
  from .macrostate import Macrostate

  strandlist = [strand_struct[:] for strand_struct in complex.structure.to_strandlist()]
  strands = complex.strands
  
  domain_start = sum([d.length
                        for d
                        in strands[strand_num].domains[:domain_num]])
  domain_end = domain_start + strands[strand_num].domains[domain_num].length
  n = 0
  for strand_n, strand_struct in enumerate(strandlist):
    for i, bound in enumerate(strand_struct):
      is_domain = (strand_n == strand_num and domain_start <= i < domain_end)
      is_bound_to_domain = (bound != None and bound != '?'
                           and bound[0] == strand_num
                           and domain_start <= bound[1] < domain_end)
      if not is_domain and not is_bound_to_domain:
        strandlist[strand_n][i] = '?'
      else:
        n += 1
  
  ms_complex = Complex(name = complex.name + "_({0},{1})".format(strand_num, domain_num),
                       strands = strands,
                       structure = strandlist)
  return Macrostate(name = "macrostate_" + ms_complex.name,
                    type = 'loose',
                    complex = ms_complex,
                    cutoff = int(cutoff * n))
        
def count_by_domain_macrostate(complex, cutoff):
  """ Returns a Macrostate that matches a complex such that every domain
  matches the given complex's structure to within the cutoff fraction. 
  cutoff is a fractional defect allowed over each domain. """
  from .macrostate import Macrostate

  macrostates = []
  for strand_num, strand in enumerate(complex.strands):
    for domain_num in range(len(strand.domains)):
      macrostates.append(loose_domain_macrostate(complex, strand_num, domain_num, cutoff))
  return Macrostate(name = complex.name + "_loose-domains",
                    type = 'conjunction',
                    macrostates = macrostates)
                    
## RestingSet -> Macrostate functions
def restingset_count_by_complex_macrostate(restingset, cutoff):
  """ Creates a macrostate corresponding to secondary structures that
  match the binding of one of the complexes in the given resting set, within the given cutoff.
  cutoff is a fractional defect over the entire complex. """
  from .macrostate import Macrostate

  macrostates = [count_by_complex_macrostate(complex, cutoff) for complex in restingset.complexes]
  return Macrostate(name        = "macrostate_{}".format(restingset.name),
                    type        = "disjunction",
                    macrostates = macrostates)
def restingset_count_by_domain_macrostate(restingset, cutoff):
  #print "WARNING: Multistrand may not support macrostates that are defined as a per-domain p-approximation"
  from .macrostate import Macrostate

  macrostates = [count_by_domain_macrostate(complex, cutoff) for complex in restingset.complexes]
  return Macrostate(name        = "macrostate_{}".format(restingset.name),
                    type        = "disjunction",
                    macrostates = macrostates)

## Functions on RestingSets
def get_containing_set(restingsets, complex):
  for rs in restingsets:
    if complex in rs:
      return rs
  return None
                    
## Functions on Structures
def num_wildcards(structure):
  s = sum(structure.to_strandlist(), [])
  return len([x for x in s if x == '?'])
  

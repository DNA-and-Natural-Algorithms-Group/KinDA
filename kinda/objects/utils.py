
## Sequence utilities
                  
def random_sequence(constraint, base_probs = None):
  import random
  from constraints import base_group
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
  
## Functions for Domain objects
def split_domain(domain, pos):  
  """ Splits a domain into two subdomains. """
  from domain import Domain
  
  section1 = Domain(name = domain.name + "[:{0}]".format(pos),
                    constraints = domain.constraints[:pos])
  section2 = Domain(name = domain.name + "[{0}:]".format(pos + 1),
                    constraints = domain.constraints[pos:])
      
  domain.subdomains = [section1, section2]

def equate_domains(domains):
  """This is pretty slow and inefficient."""
  from domain import Domain
  
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
    
    # Make a new Domain with constraints from both domains
    new_domain = Domain(name = "~(" + min(domain1.name, domain2.name) + ")",
                        constraints = domain1.constraints)
    new_domain.add_constraints(domain2.constraints)
    
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
  from structure import Structure
  
  strandlist = structure.structure[:]
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
  from macrostate import Macrostate
  
  if macrostate.type == Macrostate.types['conjunction']:
    return list(set(sum([get_dependent_complexes(m) for m in macrostate.macrostates], [])))
  else:
    return [macrostate.complex]
    
## Complex -> Macrostate functions
def exact_complex_macrostate(complex):
  from macrostate import Macrostate

  return Macrostate(name = "macrostate_" + complex.name,
                    type = 'exact',
                    complex = complex)
def loose_domain_macrostate(complex, strand_num, domain_num, cutoff):
  """ Creates a loose macrostate with the given domain as the region
  of interest. The cutoff is a fractional defect within this domain. The
  results are not meaningful if the domain is not completely bound or 
  completely unbound."""
  from complex import Complex
  from structure import Structure
  from macrostate import Macrostate

  strandlist = complex.structure.to_strandlist()[:]
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
        
def similar_complex_macrostate(complex, cutoff):
  """ Returns a Macrostate that matches a complex such that every domain
  matches the given complex's structure to within the cutoff fraction. """
  from macrostate import Macrostate

  macrostates = []
  for strand_num, strand in enumerate(complex.strands):
    for domain_num in range(len(strand.domains)):
      macrostates.append(loose_domain_macrostate(complex, strand_num, domain_num, cutoff))
  return Macrostate(name = complex.name + "_loose-domains",
                    type = 'conjunction',
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
  return len(filter(lambda x: x == '?', s))
  
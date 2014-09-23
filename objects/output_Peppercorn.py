import enumerator.utils as enum
import domain, strand, complex

def to_Peppercorn_domain(domain):
  if domain.name[-1] == "*":
    name = domain.name[:-1]
    comp = True
  else:
    name = domain.name
    comp = False
  return enum.Domain(name, domain.length, comp)
def to_Peppercorn_strand(strand, domains):
  d_list = [domains[d] for d in strand.base_domains()]
  return enum.Strand(strand.name, d_list)
def to_Peppercorn_complex(complex, strands):
  s_list = [strands[strand] for strand in complex.strands]
  
  base_to_seq_dict = {}
  for strand_num, strand in enumerate(complex.strands):
    index = 0
    for domain_num, domain in enumerate(strand.base_domains()):
      for i in range(domain.length):
        base_to_seq_dict[(strand_num, index + i)] = (strand_num, domain_num)
      index += domain.length
  base_to_seq_dict[None] = None
      
  structure = []
  for strand_struct, strand in zip(complex.structure.to_strandlist(), complex.strands):
    index = 0
    substruct = []
    for domain in strand.base_domains():
      substruct.append(base_to_seq_dict[strand_struct[index]])
      index += domain.length
    structure.append(substruct)
    
  return enum.Complex(complex.name, s_list, structure)
def to_Peppercorn(domains, strands = [], complexes = []):
  """Prepares the domain-level information for this system.  To be called 
  once all the sequences, structures, etc. have been added to the spec. 
  Raises an exception if Enumerator is not successfully created."""
  enum_domains = {}
  for d in set(sum([d.base_domains() for d in domains], [])):
    enum_domains[d] = to_Peppercorn_domain(d)
  
  enum_strands = {}
  for strand in set(strands):
    enum_strands[strand] = to_Peppercorn_strand(strand, enum_domains)
  
  enum_complexes = {}
  for complex in set(complexes):
    enum_complexes[complex] = to_Peppercorn_complex(complex, enum_strands)
  
  # Print everything for debugging purposes
  print "domains = {}"
  for old_dom, dom in enum_domains.items():
    print "domains[\"%s\"] = utils.Domain(\"%s\", %d, %s)" % \
        (old_dom.name, dom.identity, dom.length, dom.is_complement)
  
  print "strands = {}"
  for old_strand, strand in enum_strands.items():
    print "strands[\"%s\"] = utils.Strand(\"%s\", %s)" % \
        (old_strand.name, strand.name, "[" + ', '.join(["domains[\"%s\"]" % d.name for d in strand.domains]) + "]")
  
  print "complexes = {}"
  for old_complex, c in enum_complexes.items():
    print "complexes[%s] = utils.Complex(\"%s\", %s, %s)" % \
        (c.name, c.name, "[" + ', '.join(["strands[\"%s\"]" % s.name for s in c.strands]) + "]", c.structure)
  
  return (enum_domains, enum_strands, enum_complexes)
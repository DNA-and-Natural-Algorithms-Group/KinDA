from imports import peppercornhome

### TODO: Add conversion functionality to/from RestingSets

######################################
## Conversion to Peppercorn objects ##
######################################
#
# The following functions may be used to convert a list of domains, strands,
# complexes, and reactions to their Peppercorn equivalents, and may be used as
# input to Peppercorn directly with
#  enumerator.Enumerator(domains, strands, complexes)
# where domains, strands, and complexes are lists of the Peppercorn objects
# created by this system. Note that in the conversion, any composite domains
# (e.g. PIL-style supersequences) are lost because Peppercorn does not support
# composite domain definitions.
#
#
def to_Peppercorn_domain(domain):
  """ Converts a DNAObjects.Domain object to the equivalent Peppercorn
  construct. If the domain name ends in an asterisk (*) it is assumed
  to be a complementary domain. Otherwise, it is assumed to be
  non-complementary. """  
  import enumerator.utils as enum
  
  # Determine if complementary
  if domain.name[-1] == "*":
    name = domain.name[:-1]
    comp = True
    constraints = domain.constraints.complement # Input constraints for non-complement
  else:
    name = domain.name
    comp = False
    constraints = domain.constraints
  # Return Peppercorn domain
  peppercorn_domain = enum.Domain(name, domain.length, comp, constraints)
  return peppercorn_domain
  
def to_Peppercorn_strand(strand, domains):
  """ Converts a DNAObjects.Strand object to the equivalent Peppercorn
  construct. This function requires as an argument a dict mapping
  DNAObject.Domain objects to their equivalent Peppercorn Domains. """
  import enumerator.utils as enum
  
  # Make list of Peppercorn Domains from the original domains list
  d_list = [domains[d] for d in strand.base_domains()]
  # Return Peppercorn Strand
  return enum.Strand(strand.name, d_list)
  
def to_Peppercorn_complex(complex, strands):
  """ Converts a DNAObjects.Complex object to the equivalent Peppercorn
  construct. This function requires a dict mapping DNAObject Strand objects
  to their equivalent Peppercorn Strands. """
  import enumerator.utils as enum
  
  # Make list of Peppercorn Strands for this complex
  s_list = [strands[strand] for strand in complex.strands]
  
  ## Determine domain-level strandlist structure representation
  # Map (strand_num, index) tuples to (strand_num, domain_num)
  base_to_seq_dict = {}
  for strand_num, strand in enumerate(complex.strands):
    index = 0
    for domain_num, domain in enumerate(strand.base_domains()):
      for i in range(domain.length):
        base_to_seq_dict[(strand_num, index + i)] = (strand_num, domain_num)
      index += domain.length
  base_to_seq_dict[None] = None
      
  # Create new strandlist representation
  structure = []
  for strand_struct, strand in zip(complex.structure.to_strandlist(), complex.strands):
    index = 0
    substruct = []
    for domain in strand.base_domains():
      substruct.append(base_to_seq_dict[strand_struct[index]])
      index += domain.length
    structure.append(substruct)
    
  # Return Peppercorn Complex
  return enum.Complex(complex.name, s_list, structure)
  
def to_Peppercorn_reaction(reaction, complexes):
  """ Converts a DNAObjects.Reaction object to the equivalent Peppercorn
  construct. This function requires a dict mapping DNAObject Complex objects
  to their equivalent Peppercorn Complexes. """
  import enumerator.utils as enum
  
  # Compute ReactionPathway field values
  name = reaction.name
  reactants = list(reaction.reactants)
  products = list(reaction.products)
  
  # Return new ReactionPathway object
  return enum.ReactionPathway(name, reactants, products)
  
def to_Peppercorn(*args, **kargs):
  """Converts DNAObjects objects to Peppercorn objects.
  The following objects will be recognized in the key list and converted:
    - domains
    - strands
    - complexes
    - reactions
  Other supplied objects are ignored.
  The output is in the form of a dictionary mapping the keys
  'domains', 'strands', 'complexes', and 'reactions' to a list of tuples
  of the form (object, converted_object). This may be iterated through
  directly or converted into a dict for object lookup.
  """
  ## Extract all objects to be converted
  reactions = set(kargs.get('reactions', []))
  complexes = set(sum([list(r.reactants | r.products) for r in reactions], [])
                  + kargs.get('complexes', []))
  strands = set(sum([c.strands for c in complexes], [])
                + kargs.get('strands', []))
  domains = set(sum([s.base_domains() for s in strands], [])
                + sum([d.base_domains() for d in kargs.get('domains', [])], []))
                
  ## Convert domains, strands, complexes, and reactions
  enum_domains = {}
  for d in set(sum([d.base_domains() for d in domains], [])):
    enum_domains[d] = to_Peppercorn_domain(d)
  
  enum_strands = {}
  for strand in set(strands):
    enum_strands[strand] = to_Peppercorn_strand(strand, enum_domains)
  
  enum_complexes = {}
  for complex in set(complexes):
    enum_complexes[complex] = to_Peppercorn_complex(complex, enum_strands)
    
  enum_reactions = {}
  for reaction in set(reactions):
    enum_reactions[reaction] = to_Peppercorn_reaction(reaction, enum_complexes)
  
  ## Return in a dict
  results = dict()
  results['domains'] = enum_domains.items()
  results['strands'] = enum_strands.items()
  results['complexes'] = enum_complexes.items()
  results['reactions'] = enum_reactions.items()
  return results
#
###############################


  
########################################
## Conversion from Peppercorn objects ##
########################################
#
def from_Peppercorn_domain(domain):
  """ Converts a Peppercorn Domain object to an equivalent
  DNAObjects Domain object, with the given constraints. If no constraints
  are supplied, the 'N' constraint is applied to all bases in the domain. """
  import dnaobjects as dna
  
  if domain.sequence == None: constraints = dna.Constraints("N" * domain.length)
  else: constraints = dna.Constraints(domain.sequence)
  print domain, domain.sequence, constraints
  
  if domain.is_complement:
    return dna.Domain(name = domain.identity, constraints = constraints.complement).complement
  else:
    return dna.Domain(name = domain.name, constraints = constraints)
  
def from_Peppercorn_strand(strand, domains):
  """ Converts a Peppercorn Strand object to an equivalent DNAObjects Strand.
  This function requires as an argument a dict mapping Peppercorn Domain
  objects to their equivalent DNAObject Domains. """
  import dnaobjects as dna
  
  d_list = [domains[d] for d in strand.domains]
  return dna.Strand(name = strand.name, domains = d_list)
  
def from_Peppercorn_complex(complex, strands):
  """ Converts a Peppercorn Complex object to an equivalent DNAObjects Complex.
  This function requires as an argument a dict mapping Peppercorn Strand
  objects to their equivalent DNAObject Strands. """
  import dnaobjects as dna
  
  s_list = [strands[s] for s in complex.strands]
  return dna.Complex(name = complex.name,
                     strands = s_list,
                     structure = complex.structure)
                     
def from_Peppercorn_reaction(reaction, complexes):
  """ Converts a Peppercorn ReactionPathway object to an equivalent DNAObjects
  Reaction object. This function requires as an argument a dict mapping
  Peppercorn Complexes to their equivalent DNAObject Complexes. """
  import dnaobjects as dna
  
  reactants = [complexes[r] for r in reaction.reactants]
  products = [complexes[p] for p in reaction.products]
  return dna.Reaction(name = reaction.name,
                      reactants = reactants,
                      products = products)
                      
def from_Peppercorn_restingset(restingset, complexes):
  """ Converts a Peppercorn RestingSet object to an equivalent DNAObjects
  RestingSet object. This function requires as an argument a dict mapping
  Peppercorn Complexes to their equivalent DNAObjects complexes. """
  import dnaobjects as dna
  
  rs_complexes = [complexes[c] for c in restingset.complexes]
  return dna.RestingSet(name = restingset.name, complexes = rs_complexes)
                     
def from_Peppercorn(*args, **kargs):
  """ Converts a system of Peppercorn objects to the equivalent DNAObject
  counterparts. The following objects are recognized in the key arguments
  list and will be converted:
    - domains
    - strands
    - complexes
    - reactions
    - restingsets
  Other arguments are ignored. """
  import enumerator.utils as enum
  
  ## Extract all objects to be converted
  restingsets = set(kargs.get('restingsets', []))
  reactions = set(kargs.get('reactions', []))
  complexes = set(sum([r.reactants + r.products for r in reactions], [])
                  + sum([rs.complexes for rs in restingsets], [])
                  + kargs.get('complexes', []))
  strands = set(sum([c.strands for c in complexes], [])
                + kargs.get('strands', []))
  domains = set(sum([s.domains for s in strands], [])
                + kargs.get('domains', []))
                
  ## Convert objects
  dna_domains = {}
  for d in domains:
    dna_domains[d] = from_Peppercorn_domain(d)
  
  dna_strands = {}
  for s in strands:
    dna_strands[s] = from_Peppercorn_strand(s, dna_domains)
    
  dna_complexes = {}
  for c in complexes:
    dna_complexes[c] = from_Peppercorn_complex(c, dna_strands)
    
  dna_reactions = {}
  for r in reactions:
    dna_reactions[r] = from_Peppercorn_reaction(r, dna_complexes)
    
  dna_restingsets = {}
  for rs in restingsets:
    dna_restingsets[rs] = from_Peppercorn_restingset(rs, dna_complexes)
    
  ## Make return value
  result = dict()
  result['domains'] = dna_domains.items()
  result['strands'] = dna_strands.items()
  result['complexes'] = dna_complexes.items()
  result['reactions'] = dna_reactions.items()
  result['restingsets'] = dna_restingsets.items()
  return result

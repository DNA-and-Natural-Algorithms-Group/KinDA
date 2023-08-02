import itertools as it

import peppercornenumerator.objects as enumobj

from . import dnaobjects as dna

### Possible future TODOs:
###   Add conversion functionality to Peppercorn resting sets

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


def to_Peppercorn_domain(domain):
  """
  Converts a DNAObjects.Domain object to the equivalent Peppercorn construct. If
  the domain name ends in an asterisk (*) it is assumed to be a complementary
  domain. Note that the sequence of the domain is lost in the conversion, as
  Peppercorn does not store this information.
  """
  # Return Peppercorn domain
  peppercorn_domain = enumobj.PepperDomain(domain.name, length=domain.length)
  return peppercorn_domain
  

def to_Peppercorn_complex(complex, domains):
  """
  Converts a DNAObjects.Complex object to the equivalent Peppercorn construct.
  This function requires a dict mapping DNAObject Strand objects to their
  equivalent Peppercorn Strands.
  """
  nl_structure = complex.structure.to_dotparen()

  # Compute PepperComplax domain sequence and structure
  # The domain sequence is a list of domains, interpersed with '+' elements indicating nicks
  # The structure is a list of '(', '.', ')', and '+', with '+' indicating nicks
  enum_domains = []
  dl_structure = ""
  idx = 0
  for d_list in complex.base_domains():
    for d in d_list:
      enum_domains.append(domains[d])
      dl_structure += nl_structure[idx]
      idx += d.length
    enum_domains.append('+')
    dl_structure += '+'
    idx += 1
  enum_domains = enum_domains[:-1] # remove last +
  dl_structure = list(dl_structure[:-1]) # remove last +

  # Create and return equivalent PepperComplex object
  return enumobj.PepperComplex(enum_domains, dl_structure, name = complex.name)
  
 
def to_Peppercorn_reaction(reaction, complexes):
  """
  Converts a DNAObjects.Reaction object to the equivalent Peppercorn construct.
  This function requires a dict mapping DNAObject Complex objects to their
  equivalent Peppercorn Complexes.
  """
  # Compute PepperReaction field values
  name = reaction.name
  reactants = [complexes[c] for c in reaction.reactants]
  products = [complexes[c] for c in reaction.products]
  
  # Return new ReactionPathway object
  return enumobj.PepperReaction(name, reactants, products)
  

def to_Peppercorn(*args, **kargs):
  """
  Converts DNAObjects objects to Peppercorn objects.

  The following objects will be recognized in the key list and converted:
    - domains
    - complexes
    - reactions
  Other supplied objects are ignored.

  The output is in the form of a dictionary mapping the keys 'domains',
  'complexes', and 'reactions' to a list of tuples of the form (object,
  converted_object). This may be iterated through directly or converted into a
  dict for object lookup.
  """
  ## Clear previously created Peppercorn objects from memory
  ## to prevent duplication errors.
  enumobj.clear_memory()
 
  ## Extract all objects to be converted
  reactions = set(kargs.get('reactions', []))
  complexes = set(sum([list(r.reactants + r.products) for r in reactions], [])
                  + kargs.get('complexes', []))
  domains = set([d for c in complexes for d_list in c.base_domains() for d in d_list]
                + sum([d.base_domains() for d in kargs.get('domains', [])], []))
                
  ## Convert domains, strands, complexes, and reactions
  enum_domains = {}
  for d in set(sum([d.base_domains() for d in domains], [])):
    enum_domains[d] = to_Peppercorn_domain(d)
  
  enum_complexes = {}
  for complex in set(complexes):
    enum_complexes[complex] = to_Peppercorn_complex(complex, enum_domains)
    
  enum_reactions = {}
  for reaction in set(reactions):
    enum_reactions[reaction] = to_Peppercorn_reaction(reaction, enum_complexes)
  
  ## Return in a dict
  results = dict()
  results['domains'] = list(enum_domains.items())
  results['complexes'] = list(enum_complexes.items())
  results['reactions'] = list(enum_reactions.items())
  return results


########################################
## Conversion from Peppercorn objects ##
########################################


def from_Peppercorn_domain(domain, seq = None):
  """
  Converts a Peppercorn PepperDomain object to an equivalent DNAObjects Domain
  object, with the given sequence. If no sequence is supplied, the 'N'
  constraint is applied to all bases in the domain.
  """
  if seq is None:  sequence = dna.Sequence("N"* domain.length)
  else:  sequence = dna.Sequence(seq)

  if domain.is_complement:
    return dna.Domain(name = domain.cname, sequence = sequence.complement).complement
  else:
    return dna.Domain(name = domain.name, sequence = sequence)


def from_Peppercorn_strand(strand, domains, name = None):
  """
  Converts a Peppercorn 'strand' (i.e. a tuple of PepperDomain objects) into an
  equivalent DNAObjects Strand object. Note that because Peppercorn no longer
  implements an object representing strands, this tuple representation is the
  closest we have. This means that strand names are not preserved when
  converting to Peppercorn and back.
  """
  if name is None:
    return dna.Strand(domains = [domains[d] for d in strand])
  else:
    return dna.Strand(name = name, domains = [domains[d] for d in strand])


def from_Peppercorn_complex(complex, strands):
  """
  Converts a Peppercorn PepperComplex object to an equivalent DNAObjects
  Complex. This function requires as an argument a dict mapping tuples of
  Peppercorn PepperDomains to their equivalent DNAObject Strand objects. Because
  Peppercorn does not currently implement a Strand object, the PepperDomain
  tuples give us a consistent representation between PepperComplexes. It is
  important that, after the conversion, different DNAObject Complexes share
  Strand objects, or Multistrand simulations may not work.
  """
  pepper_s_list = [tuple(y) for x,y in
                   it.groupby(complex.sequence, lambda char: char=='+') if not x]
  s_list = [strands[pepper_s] for pepper_s in pepper_s_list]
  return dna.Complex(name = complex.name,
                     strands = s_list,
                     structure = ''.join(complex.structure))
                     

def from_Peppercorn_reaction(reaction, complexes):
  """
  Converts a Peppercorn PepperReaction object to an equivalent DNAObjects
  Reaction object. This function requires as an argument a dict mapping
  Peppercorn PepperComplex objects to their equivalent DNAObject Complex. Note:
  Peppercorn PepperReaction objects do not store reaction names, so the new
  DNAObjects Reaction will have an automatically generated name.
  """
  reactants = [complexes[r] for r in reaction.reactants]
  products = [complexes[p] for p in reaction.products]
  return dna.Reaction(reactants = reactants,
                      products = products)
                      

def from_Peppercorn_restingset(restingset, complexes):
  """
  Converts a Peppercorn RestingSet object to an equivalent DNAObjects RestingSet
  object. This function requires as an argument a dict mapping Peppercorn
  Complexes to their equivalent DNAObjects complexes.
  """
  rs_complexes = [complexes[c] for c in restingset.complexes]
  return dna.RestingSet(name = restingset.name, complexes = rs_complexes)
                     

def from_Peppercorn_restingset_reaction(reaction, restingsets):
  """
  Converts a condensed Peppercorn PepperReaction object to an equivalent
  DNAObjects RestingSetReaction object. This function requires as an argument a
  dict mapping Peppercorn RestingSet objects to equivalent DNAObject
  RestingSets.
  """
  if reaction.rtype != 'condensed':
    print("KinDA: WARNING: Attempted to convert non-condensed Peppercorn "
          "reaction into a RestingSetReaction")
  
  reactants = [restingsets[r] for r in reaction.reactants]
  products = [restingsets[p] for p in reaction.products]
  return dna.RestingSetReaction(reactants = reactants, products=products)


def from_Peppercorn(*args, **kargs):
  """
  Converts a system of Peppercorn objects to the equivalent DNAObject
  counterparts. The following objects are recognized in the key arguments list
  and will be converted:
    - domains
    - complexes
    - reactions
    - restingsets
    - restingsetreactions
  Other arguments are ignored. Note that Peppercorn does not currently have a
  representation for strands, so we do not handle them here. Because the
  conversion to Peppercorn objects is lossy (in particular, domains lose
  sequence information and strands lose name information), we allow two
  additional parameters:
    - domain_info: a dict mapping domain names to dict item lists of
                   additional domain info
    - strand_info: a dict mapping tuples of domain names to dict item lists of
                   strand info
  """
  ## Extract all objects to be converted
  restingsetreactions = set(kargs.get('restingsetreactions', []))
  restingsets = set(it.chain(*it.chain(
    [kargs.get('restingsets', [])],
    *([r.reactants, r.products] for r in restingsetreactions))))
  reactions = set(kargs.get('reactions', []))
  complexes = set(it.chain(*it.chain(
    *([r.reactants, r.products] for r in reactions),
    (rs.complexes for rs in restingsets),
    [kargs.get('complexes', [])])))
  strands = set([tuple(strand) for c in complexes for nick,strand in it.groupby(c.sequence, lambda v:v=='+') if not nick])
  domains = set([d for strand in strands for d in strand]
                + kargs.get('domains', []))
  domain_info = kargs.get('domain_info', {})
  strand_info = kargs.get('strand_info', {})
                
  ## Convert objects
  dna_domains = {}
  for d in domains:
    d_info = dict(domain_info.get(d.name, []))
    dna_domains[d] = from_Peppercorn_domain(d, **d_info)

  dna_strands = {}
  for s in strands:
    s_key = tuple([dna_domains[d].name for d in s])
    s_info = dict(strand_info.get(s_key, []))
    dna_strands[s] = from_Peppercorn_strand(s, dna_domains, **s_info)

  dna_complexes = {}
  for c in complexes:
    dna_complexes[c] = from_Peppercorn_complex(c, dna_strands)
    
  dna_reactions = {}
  for r in reactions:
    dna_reactions[r] = from_Peppercorn_reaction(r, dna_complexes)
    
  dna_restingsets = {}
  for rs in restingsets:
    dna_restingsets[rs] = from_Peppercorn_restingset(rs, dna_complexes)

  dna_restingsetreactions = {}
  for rs_rxn in restingsetreactions:
    dna_restingsetreactions[rs_rxn] = from_Peppercorn_restingset_reaction(
      rs_rxn, dna_restingsets)
    
  ## Make return value
  result = dict()
  result['domains'] = list(dna_domains.items())
  result['complexes'] = list(dna_complexes.items())
  result['reactions'] = list(dna_reactions.items())
  result['restingsets'] = list(dna_restingsets.items())
  result['restingsetreactions'] = list(dna_restingsetreactions.items())
  return result

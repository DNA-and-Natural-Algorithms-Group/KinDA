## IMPORTS
import re

import dnaobjects as dna

## GLOBALS
# Helpful regular expressions
space   = r"\s*" # 0 or more spaces or tabs
word    = r"[\w\-*]+"
words_s = r"(?:" + word + space + r")+"             # space-separated words
words_p = r"(?:\+" + space + word + space + r")*"   # plus-separated words
seq     = r"[ACGTSWMKBVDHRYN]+"
dec     = r"(?:(?:(?:\d+(?:\.\d*)?)|\.\d+)(?:[eE][+-]?[0-9]+)?|inf)"
c_unit  = r"[nmup]?M"
ic_unit  = r"(?:/[nmup]?M)"
it_unit  = r"(?:/s(?!M)|/sec(?!M)|/min(?!M)|/m(?!M))"
unit    = "(?:" + ic_unit + "?" + space + it_unit + ")"
details = r"\[.*\]"
cmax    = r"(?:cmax" + space + "=" + space + dec + space + c_unit + r")"
opt     = r"\d+nt"
kin     = r"\[kinetic\]"
dummy   = r"(?:\[dummy\])?"
d_paren = r"[().+]+"
num     = r"\d+"

## FUNCTIONS
def read_PIL(filename):
  """Parses the specified PIL file and returns the circuit design objects."""

  spec_file = open(filename)
  domains = {}
  strands = {}
  complexes = []
  # Parse each line with regular expressions
  for l in spec_file:
    #print l,
    line = strip_comment(l).strip()
    if line == "":
      continue
      
    directive = line.split()[0]
    
    # If sequence declaration
    if directive == "sequence":
      d = parseSequenceDirective(line)
      domains[d.name] = d
      domains[d.complement.name] = d.complement
    elif directive == "sup-sequence":
      d = parseSupseqDirective(line, domains)
      domains[d.name] = d
      domains[d.complement.name] = d.complement
    elif directive == "strand":
      s = parseStrandDirective(line, domains)
      strands[s.name] = s
      strands[s.complement.name] = s.complement
    elif directive == "structure":
      complexes.append(parseStructDirective(line, strands))
    elif directive == "kinetic": pass
      #parseKineticDirective(line, spec)
    elif directive == "equal":
      parseEqualDirective(line, domains)
    elif directive == "noninteracting": pass
      #parseNoninteractingDirective(line, spec)
    else:
      # Unknown directive on this line; print error message and continue
      print >> sys.stderr, "Warning: Bad syntax on this line:\n%s\nContinuing anyway.\n" % l
      
  spec_file.close()
  
  return (domains.values(), strands.values(), complexes)
  

def parseSequenceDirective(line):
  # Syntax:
  #   sequence <name> = <nucleotides> : <seq_length>
  regex = re.compile("sequence" + space + "(" + word + ")" + space + "="
    + space + "(" + seq + ")" + space + ":" + space + num)
  match_obj = regex.match(line)
  if match_obj != None:
    name = match_obj.group(1)
    nucleotides = match_obj.group(2)
    # print name, "=", nucleotides ####
    
    # Create/add the new Sequence to the Specification
    domain = dna.Domain(name = name, constraints = nucleotides)
    return domain
  else:
    print >> sys.stderr, "Invalid sequence directive:\n%s" % line
    
    
def parseSupseqDirective(line, domains):
  # Syntax:
  #   sup-sequence <name> = <seq1> [<seq2> [<seq3> ...]] : <seq_length>
  regex = re.compile("sup-sequence" + space + "(" + word + ")" + space + "=" +
    space + "(" + words_s + ")" + ":" + space + num)
  match_obj = regex.match(line)
  if match_obj != None:
    name = match_obj.group(1)
    seq_list = [domains[n] for n in match_obj.group(2).split()]
    # print name, "=", match_obj.group(2).split() ####
      
    # Create/add the new SuperSequence and its complement
    domain = dna.Domain(name = name, subdomains = seq_list)
    return domain
  else:
    print >> sys.stderr, "Invalid sup-sequence directive:\n%s" % line
    
        
def parseStrandDirective(line, domains):
  # Syntax:
  #   strand [\[dummy\]] <name> = <seq1> [<seq2> [<seq3> ...]] : <strand_length>
  regex = re.compile("strand" + space + dummy + space + "(" + word + ")" + space + "=" +
    space + "(" + words_s + ")" + ":" + space + num)
  match_obj = regex.match(line)
  if match_obj != None:
    name = match_obj.group(1)
    seq_list = [domains[n] for n in match_obj.group(2).split()]
    # print name, "=", match_obj.group(2).split() ####
    
    strand = dna.Strand(name = name, domains = seq_list)
    return strand
  else:
    print >> sys.stderr, "Invalid strand directive:\n%s" % line
    
      
def parseStructDirective(line, strands):
  # Syntax:
  #   structure [\[<options>\]] <name> = <strand1> [+ <strand2> [+ <strand3> ...]] : <dot_parens>
  """def parseOptions(opt_str):
    if opt_str == None:
      return (None, None)
  
    opt_list = [opt.strip() for opt in opt_str.split(",")]
    re_nt = re.compile(r"(\d)nt")
    re_cmax = re.compile(r"cmax\s*=\s*(\d+.?\d*)\s*(\w+)")
    opt_value = None
    c_max = None
    for opt in opt_list:
      mo_nt = re_nt.match(opt)
      mo_cmax = re_cmax.match(opt)
      if mo_nt != None:
        opt_value = mo_nt.group(1)
      elif mo_cmax != None:
        cmax_value = float(mo_cmax.group(1))
        cmax_unit = mo_cmax.group(2)
        c_max = convert_to_M(cmax_value, cmax_unit)
      else:
        print >> sys.stderr, "Invalid structure option:%s. Ignoring..." % opt
    return (opt_value, c_max)"""
    
  regex = re.compile("structure" + space + "(" + details + ")?" + space + "(" + word + ")" + space + "=" +
    space + "(" + word + space + words_p + ")" + ":" + space + "(" + d_paren + ")")
  match_obj = regex.match(line)
  if match_obj != None:
    #opt_value, c_max = parseOptions(match_obj.group(1)[1:-1])
    name = match_obj.group(2)
    strand_list = [strands[strand_name.strip()] for strand_name in match_obj.group(3).split("+")]
    binding = match_obj.group(4)
    # print "opt_val=%s, cmax=%s, name=%s, strandlist=%s, binding=%s" % (opt_value, c_max, name, strand_list, binding) ####
    
    # Create/add the new Structure
    complex = dna.Complex(name = name, strands = strand_list, structure = binding) # we're ignoring c_max, opt_value which is not okay
    return complex
  else:
    print >> sys.stderr, "Invalid structure directive:\n%s" % line
    
      
def parseKineticDirective(line):
  # Syntax:
  #   kinetic [\[<num> <units> < k < <num> <units>\]] <struct1> [+ <struct2> ...] -> <struct3> [+ <struct4> ...]
  k_range = (r"(?:\[" + space + "(" + dec + ")" + space + "(" + unit + ")" + space + "<" +
    space + "k" + space + "<" + space + "(" + dec + ")" + space + "(" + unit + ")" + space + r"\])")
  regex = re.compile("kinetic" + space + "(" + k_range + ")?" + space + "(" + word + space + words_p + ")->" +
    space + "(" + word + space + words_p + ")")
  match_obj = regex.match(line)
  if match_obj != None:
    if match_obj.group(1) != None:
      lower_k = float(match_obj.group(2))
      lower_units = ["/"+u for u in match_obj.group(3).split("/")]
      upper_k = float(match_obj.group(4))
      upper_units = ["/"+u for u in match_obj.group(5).split("/")]
      
      lower_cunit = lower_tunit = None
      for u in lower_units:
        if re.match(r"(/[nmup]?M)", u):
          lower_cunit = u
        if re.match(r"(/s(?!M)|/sec(?!M)|/min(?!M)|/m(?!M))", u):
          lower_tunit = u
      
      upper_cunit = upper_tunit = None
      for u in upper_units:
        if re.match(r"(/[nmup]?M)", u):
          upper_cunit = u
        if re.match(r"(/s(?!M)|/sec(?!M)|/min(?!M)|/m(?!M))", u):
          upper_tunit = u
      
      lower_k = convert_to_Ms(lower_k, lower_cunit, lower_tunit)
      if lower_k == float('inf'):
        lower_k = INFINITY
      upper_k = convert_to_Ms(upper_k, upper_cunit, upper_tunit)
      if upper_k == float('inf'):
        upper_k = INFINITY
      krange = [lower_k, upper_k]
    else:
      krange = [None, None]
      
    ins = [spec.structure_dict[s.strip()] for s in match_obj.group(6).split("+")]
    outs = [spec.structure_dict[s.strip()] for s in match_obj.group(7).split("+")]
    
    # print "krange = %s, ins = %s, outs = %s" % (krange, match_obj.group(6).split("+"), match_obj.group(7).split("+")) ####
      
    #spec.add_kinetic(Kinetic(ins, outs, *krange))
    # DON'T RETURN ANYTHING BECAUSE WE DON'T HANDLE THESE STATEMENTS YET :(
  else:
    print >> sys.stderr, "Invalid kinetic directive:\n%s" % line
    
    
def parseEqualDirective(line, domains):
  # Syntax:
  #   equal <seq1> [<seq2> [<seq3> ...]]
  regex = re.compile("equal" + space + "(" + words_s + ")")
  match_obj = regex.match(line)
  if match_obj != None:
    seq_list = set([domains[name] for name in match_obj.group(1).split()])
    dna.utils.equate_domains(seq_list)
  else:
    print >> sys.stderr, "Invalid equal directive:\n%s" % line
      
      
def parseNoninteractingDirective(line):
  # Syntax:
  #   noninteracting \[kinetic\] <struct1> [<struct2> ...]
  regex = re.compile("noninteracting" + space + kin + space + "(" + words_s + ")")
  match_obj = regex.match(line)
  if match_obj != None: pass
    # print match_obj.groups() ####
    #spec.noninteracting = [spec.structure_dict[name] for name in match_obj.group(1).split()] # does this work?
    # DON'T RETURN ANYTHING BECAUSE WE DON'T HANDLE THESE STATEMENTS YET :(
  else:
    print >> sys.stderr, "Invalid noninteracting directive:\n%s" % line
  

def strip_comment(line):
  try:
    return line[:line.index("#")]
  except ValueError:
    return line
    
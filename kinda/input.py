
from __future__ import print_function
import itertools as it

from dsdobjects.parser import parse_kernel_string, parse_kernel_file
from dsdobjects.parser import parse_pil_string, parse_pil_file
from dsdobjects.parser import ParseException, PilFormatError

from .objects import dnaobjects as dna

def resolve_loops(loop):
    """ Return a sequence, structure pair from kernel format with parenthesis. """
    sequen = []
    struct = []
    for dom in loop :
        if isinstance(dom, str):
            sequen.append(dom)
            if dom == '+' :
                struct.append('+')
            else :
                struct.append('.')
        elif isinstance(dom, list):
            struct[-1] = '('
            old = sequen[-1]
            se, ss = resolve_loops(dom)
            sequen.extend(se)
            struct.extend(ss)
            sequen.append(old + '*' if old[-1] != '*' else old[:-1])
            struct.append(')')
    return sequen, struct

def read_pil(data, is_file = False, composite = False):
    """ Read PIL file format.

    Use dsdobjects parser to extract information. Load kinda.objects.

    Args:
        data (str): Is either the PIL file in string format or the path to a file.
        is_file (bool): True if data is a path to a file, False otherwise
    """
    if is_file :
        parsed_file = parse_pil_file(data)
    else :
        parsed_file = parse_pil_string(data)

    domains = {'+' : '+'} # saves some code
    strands = {}
    get_strand = {}

    complexes = {}
    resting = {}
    con_reactions = []
    det_reactions = []
    for line in parsed_file :
        name = line[1]
        if line[0] == 'dl-domain':
            raise PilFormatError('KinDA needs nucleotide level information.')

        elif line[0] == 'sl-domain':
            if len(line) == 4:
                if int(line[3]) != len(line[2]):
                    raise PilFormatError("Sequence/Length information inconsistent {} vs ().".format(
                        line[3], len(line[2])))

            sequence = dna.Sequence(line[2])

            if name[-1] == '*':
                # This will be possible, sooner or later. But we have to make sure the 
                # kinda.objects can handle it.
                raise NotImplementedError
            else :
                dom = dna.Domain(name = name, sequence = line[2])
                domains[dom.name] = dom
                domains[dom.complement.name] = dom.complement

        elif line[0] == 'composite-domain':
            domain_list = map(lambda x: domains[x], line[2])

            d = dna.Domain(name = name, subdomains = domain_list)
            domains[d.name] = d
            domains[d.complement.name] = d.complement

            # if it is a strand...
            s = dna.Strand(name = name, domains = domain_list)
            strands[s.name] = s
            strands[s.complement.name] = s.complement
            
            get_strand[tuple(domain_list)] = s


        elif line[0] == 'strand-complex':
            strand_list = map(lambda x: strands[x], line[2])
            structure = line[3].replace(' ','')
            cplx = dna.Complex(name = name, strands = strand_list, structure = structure ) 
            complexes[cplx.name] = cplx
 
        elif line[0] == 'kernel-complex':
            sequence, structure = resolve_loops(line[2])

            # Replace names with domain objects.
            try :
                sequence = map(lambda d : domains[d], sequence)
            except KeyError:
                raise PilFormatError("Cannot find domain: {}.".format(d))
            
            current = []
            strand_list = []
            for d in sequence + ['+']:
                if isinstance(d, dna.Domain):
                    current.append(d)
                else:
                    if tuple(current) not in get_strand:
                        sname = '_'.join(map(str, current))
                        s = dna.Strand(name = sname, domains = current)
                        strands[s.name] = s
                        strands[s.complement.name] = s.complement
                        get_strand[tuple(current)] = s

                    strand_list.append(get_strand[tuple(current)])
                    current = []

            cplx = dna.Complex(name = name, strands = strand_list, structure = ''.join(structure))
            complexes[cplx.name] = cplx

        elif line[0] == 'resting-macrostate':
            cplxs = map(lambda c : complexes[c], line[2])
            resting[name] = dna.RestingSet(name = name, complexes = cplxs)

        elif line[0] == 'reaction':
            rtype = line[1][0][0] if line[1] != [] and line[1][0] != [] else None

            assert rtype is not None
            
            if rtype == 'condensed' :
                reactants = map(lambda c : resting[c], line[2])
                products  = map(lambda c : resting[c], line[3])
                con_reactions.append(
                        dna.RestingSetReaction(reactants = reactants, products = products))
            else :
                reactants = map(lambda c : complexes[c], line[2])
                products  = map(lambda c : complexes[c], line[3])
                det_reactions.append(
                        dna.Reaction(reactants = reactants, products = products))

        else :
            print('# Ignoring keyword: {}'.format(line[0]))

    # Make sure the reverse reaction between every pair of reactants is
    # included. These unproductive reactions will be important stop states for
    # Mulstistrand simulations. 
    reactant_pairs = it.product(resting.values(), resting.values())
    for reactants in reactant_pairs:
        con_reactions.append(dna.RestingSetReaction(reactants = reactants, 
                                                    products  = reactants))

    return complexes.values(), det_reactions, resting.values(), con_reactions


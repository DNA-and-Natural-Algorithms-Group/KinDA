"""
Provides utility functions for exporting/importing data from a KinDA session.
"""

import json
import pickle
import itertools as it

import numpy as np

from dsdobjects.dsdparser import parse_pil_string, parse_pil_file
from peppercornenumerator.input import PilFormatError

from .. import __version__, System
from .. import objects as dna
from . import io_PIL


###############################
#       Export Utilities      #
###############################

def export_data(sstats, filepath, use_pickle = False):
  """ Exports data of this KinDA object so that it can be imported in a later Python session.
  Does not export the entire KinDA object (only the XXXStats data that has been collected).
  Data is exported in JSON format.
  The following constructs are exported:
    - domains
    - strands
    - complexes
    - resting sets
    - reactions
    - resting-set reactions
    - resting-set statistics objects
    - resting-set-reaction statistics objects
  Implemented, but badly and fragily dependent on KinDA object implementation """
  ## Extract all objects to be exported
  rs_reactions = sstats._condensed_reactions
  restingsets = sstats._restingsets
  reactions = sstats._detailed_reactions
  complexes = sstats._complexes
  strands = set(sum([c.strands for c in complexes], []))
  domains = set(sum([s.base_domains() for s in strands], []))

  # Note: We destroy the domain hierarchy and only store "base" domains.
  domain_to_id = {dom: 'dom{0}_{1}'.format(dom.id, dom.name) for dom in domains}
  # id #s are guaranteed to be unique within a Python session
  domain_to_dict = {
      d_id: {'name': d.name, 'sequence': str(d.sequence)} for d,d_id in domain_to_id.items()}

  strand_to_id = {strand: 'strand{0}_{1}'.format(strand.id, strand.name) for strand in strands}
  strand_to_dict = {}
  for s, s_id in strand_to_id.items():
    strand_domains = [domain_to_id[d] for d in s.base_domains()]
    strand_to_dict[s_id] = {'name': s.name, 'domains': strand_domains}

  complex_to_id = {cpx: 'cpx{0}_{1}'.format(cpx.id, cpx.name) for cpx in complexes}
  complex_to_dict = {}
  for c, c_id in complex_to_id.items():
    cpx_strands = [strand_to_id[s] for s in c.strands]
    complex_to_dict[c_id] = {
        'name': c.name, 'strands': cpx_strands, 'structure': c.structure.to_dotparen()}

  rs_to_id = {rs: 'rs{0}_{1}'.format(rs.id, rs.name) for rs in restingsets}
  rs_to_dict = {}
  for rs, rs_id in rs_to_id.items():
    rs_complexes = [complex_to_id[c] for c in rs.complexes]
    rs_to_dict[rs_id] = {'name': rs.name, 'complexes': rs_complexes}

  rxn_to_id = {rxn: 'rxn{0}_{1}'.format(rxn.id, str(rxn)) for rxn in reactions}
  rxn_to_dict = {}
  for r, r_id in rxn_to_id.items():
    reactants = [complex_to_id[c] for c in r.reactants]
    products = [complex_to_id[c] for c in r.products]
    rxn_to_dict[r_id] = {'name': r.name, 'reactants': reactants, 'products': products}

  rsrxn_to_id = {rsrxn: 'rsrxn{0}_{1}'.format(rsrxn.id, str(rsrxn))
                 for rsrxn in rs_reactions}
  rsrxn_to_dict = {}
  for r, r_id in rsrxn_to_id.items():
    reactants = [rs_to_id[rs] for rs in r.reactants]
    products = [rs_to_id[rs] for rs in r.products]
    rsrxn_to_dict[r_id] = {'name': r.name, 'reactants': reactants, 'products': products}

  rsstats_to_dict = {}
  for rs in restingsets:
    stats = sstats.get_stats(rs)
    rsstats_to_dict[rs_to_id[rs]] = {
        'similarity_threshold': stats.get_similarity_threshold(), 'c_max': stats.c_max}
    for c in rs.complexes:
      c_data = list(stats.get_conformation_prob_data(c.name))
      rsstats_to_dict[rs_to_id[rs]][complex_to_id[c]] = {
        'prob': '{0} +/- {1}'.format(
          stats.get_conformation_prob(c.name, 1, max_sims=0),
          stats.get_conformation_prob_error(c.name, max_sims=0)),
        'similarity_data': c_data
      }
      assert len(c_data) == stats.sampler.get_num_sims()

  rsrxnstats_to_dict = {}
  for rsrxn in rs_reactions:
    stats = sstats.get_stats(rsrxn)
    sim_data = {key: d.tolist() for key,d in stats.get_simulation_data().items()}
    if len(rsrxn.reactants) == 2:
      rsrxnstats_to_dict[rsrxn_to_id[rsrxn]] = {
        'prob': '{0} +/- {1}'.format(stats.get_prob(max_sims = 0), stats.get_prob_error(max_sims=0)),
        'kcoll': '{0} +/- {1}'.format(stats.get_kcoll(max_sims = 0), stats.get_kcoll_error(max_sims=0)),
        'k1': '{0} +/- {1}'.format(stats.get_k1(max_sims = 0), stats.get_k1_error(max_sims=0)),
        'k2': '{0} +/- {1}'.format(stats.get_k2(max_sims = 0), stats.get_k2_error(max_sims=0)),
        'simulation_data': sim_data,
        'invalid_simulation_data': stats.get_invalid_simulation_data(),
        'tag': stats.multijob_tag
      }
    elif len(rsrxn.reactants) == 1:
      rsrxnstats_to_dict[rsrxn_to_id[rsrxn]] = {
        'prob': '{0} +/- {1}'.format(stats.get_prob(max_sims = 0), stats.get_prob_error(max_sims=0)),
        'k1': '{0} +/- {1}'.format(stats.get_k1(max_sims = 0), stats.get_k1_error(max_sims=0)),
        'k2': '{0} +/- {1}'.format(stats.get_k2(max_sims = 0), stats.get_k2_error(max_sims=0)),
        'simulation_data': sim_data,
        'invalid_simulation_data': stats.get_invalid_simulation_data(),
        'tag': stats.multijob_tag
      }

  # Prepare the overall dict object to be JSON-ed
  sstats_dict = {
    'domains': domain_to_dict,
    'strands': strand_to_dict,
    'complexes': complex_to_dict,
    'resting-sets': rs_to_dict,
    'reactions': rxn_to_dict,
    'resting-set-reactions': rsrxn_to_dict,
    'resting-set-stats': rsstats_to_dict,
    'resting-set-reaction-stats': rsrxnstats_to_dict,
    'initialization_params': sstats.initialization_params,
    'version': __version__
  }

  if use_pickle :
    pickle.dump(sstats_dict, open(filepath, "wb"))
  else :
    json.dump(sstats_dict, open(filepath, 'w'))


def format_rate_units(rate, arity, molarity, time):
    """ Rate must be given in /M and /s. """
    if time == 's':
        pass
    elif time == 'm':
        rate *= 60
    elif time == 'h':
        rate *= 3600
    else :
        raise NotImplementedError('Unknown time unit: {}'.format(time))

    if molarity == 'M':
        pass
    elif molarity == 'mM':
        if arity[0] > 1:
            factor = arity[0] - 1
            rate /= (factor * 1e3)
    elif molarity == 'uM':
        if arity[0] > 1:
            factor = arity[0] - 1
            rate /= (factor * 1e6)
    elif molarity == 'nM':
        if arity[0] > 1:
            factor = arity[0] - 1
            rate /= (factor * 1e9)
    else :
        raise NotImplementedError('Unknown concentration unit: {}'.format(molarity))

    return rate


def write_pil(KindaSystem: System, fh, spurious=False, unproductive=False,
              molarity='nM', time='s', prefix=None):
    """Write the KindaSystem object into a proper *.pil format.

    Args:
        KindaSystem (:obj:`kinda.System()`): The kinda.System() object.
        fh (filehandle): Filehandle to write to.
        spurious (bool, optional): Print information about
                                   spurious complexes and reactions.
        unproductive (bool, optional): Print information about unproductive
                                       complexes and reactions.

    NOTE: Eventually, this function should return more than just the CRN. That
        means info about domains, complexes, strands (ugh), etc.
    """
    # counts number of intermediate species.
    I = 0

    out = []
    def output_string(string):
        if fh is None:
            out.append(string)
        else :
            fh.write(string)

    output_string("# File generated by KinDA-{}\n".format(__version__))

    output_string('\n# Condensed reactions\n')
    for rxn in KindaSystem.get_reactions(spurious=spurious, unproductive=unproductive):
        stats = KindaSystem.get_stats(rxn)

        k_1 = stats.get_k1(max_sims=0)
        k_1_err = stats.get_k1_error(max_sims=0)
        k_2 = stats.get_k2(max_sims=0)
        k_2_err = stats.get_k2_error(max_sims=0)

        reactants = [x.name for x in rxn.reactants]
        products  = [x.name for x in rxn.products]
        if prefix :
            inter = prefix + str(I)
        else :
            inter = '_'.join(sorted(reactants) + ['to'] + sorted(products))

        if k_1 > 0 and k_2 > 0 :
            output_string('  reaction [k1 = {:10.4e} +/- {:10.4e} {:5s}] {} -> {}\n'.format(
                format_rate_units(k_1, (len(reactants), 1), molarity, time),
                format_rate_units(k_1_err, (len(reactants), 1), molarity, time),
                '/{}'.format(molarity)*(len(reactants)-1)+'/{}'.format(time),
                ' + '.join(reactants), inter))

            output_string('  reaction [k2 = {:10.4e} +/- {:10.4e} {:5s}] {} -> {}\n'.format(
                format_rate_units(k_2, (1,len(products)), molarity, time),
                format_rate_units(k_2_err, (1,len(products)), molarity, time),
                '/{}'.format(time),
                inter, ' + '.join(products)))
        else:
            output_string('# reaction [k1 = {:10.4e} +/- {:10.4e} {:5s}] {} -> {}\n'.format(
                format_rate_units(k_1, (len(reactants), 1), molarity, time),
                format_rate_units(k_1_err, (len(reactants), 1), molarity, time),
                '/{}'.format(molarity)*(len(reactants)-1)+'/{}'.format(time),
                ' + '.join(reactants), inter))

            output_string('# reaction [k2 = {:10.4e} +/- {:10.4e} {:5s}] {} -> {}\n'.format(
                format_rate_units(k_2, (1,len(products)), molarity, time),
                format_rate_units(k_2_err, (1,len(products)), molarity, time),
                '/{}'.format(time),
                inter, ' + '.join(products)))
        I += 1

    output_string('\n# Resting macrostate probabilities\n')
    restingsets = KindaSystem.get_restingsets(spurious=spurious)
    rms_width = max(len(rms.name) for rms in restingsets)
    for rms in restingsets:
        stats = KindaSystem.get_stats(rms)

        p = 1 - stats.get_conformation_prob(None, max_sims=0)
        p_err = stats.get_conformation_prob_error(None, max_sims=0)
        temp_dep = stats.get_temporary_depletion(max_sims=0)

        output_string(f'# {rms.name:{rms_width}s} [Prob = {p: >8.4%} '
                      f'+/- {p_err: >8.4%}; Depletion = {temp_dep: >8.4%}]\n')

    # In a table, print out temporary depletion levels for each pairwise
    # combination of resting sets.
    restingsets = KindaSystem.get_restingsets(spurious=False)
    if unproductive is None:
        output_string('\n# Temporary depletion details\n')
        for rms1 in restingsets:
            rms_stats = KindaSystem.get_stats(rms1)

            for rms2 in restingsets:
                rxns = KindaSystem.get_reactions(unproductive=True, reactants=[rms1,rms2])
                assert len(rxns) == 1
                rxn_stats = KindaSystem.get_stats(rxns[0])
                depl = rms_stats.get_temporary_depletion_due_to(rxn_stats, max_sims=0)

                output_string('# {:20s} [Depletion due to {:20s} = {: >8.4%}]\n'.format(
                    rms1.name, rms2.name, depl))


###############################
#       Import Utilities      #
###############################

def import_data(filepath, use_pickle = False):
  """ Imports a KinDA object as exported in the format specified by export_data()

  Imports:
    - domains, strands, complexes, reactions, resting-sets, resting-set reactions
    - resting-set stats:
        => similarity-threshold
        => c_max (concentration maximum to calculate temporary depletion)
        => list of simulation results
            - loaded int "nupackjob"
    - resting-set reaction stats:
        => load
  """
  if use_pickle:
    sstats_dict = pickle.load(open(filepath, "rb"))
  else:
    sstats_dict = json.load(open(filepath))

  if 'version' not in sstats_dict:
    print(f"# KinDA: WARNING: Imported data file has no version number. "
          f"Assuming KinDA {__version__}.")
  elif sstats_dict['version'] != __version__:
    print(f"# KinDA: WARNING: Attempting conversion from "
          f"KinDA {sstats_dict['version']}.")
    sstats_dict = _import_data_convert_version(sstats_dict, sstats_dict['version'])

  domains = {}
  for domain_id, data in sstats_dict['domains'].items():
    domains[domain_id] = dna.Domain(name = data['name'], sequence = data['sequence'])

  strands = {}
  for strand_id, data in sstats_dict['strands'].items():
    strand_domains = [domains[d_id] for d_id in data['domains']]
    strands[strand_id] = dna.Strand(name = data['name'], domains = strand_domains)

  complexes = {}
  for complex_id, data in sstats_dict['complexes'].items():
    cpx_strands = [strands[s_id] for s_id in data['strands']]
    complexes[complex_id] = dna.Complex(name = data['name'],
        strands = cpx_strands, structure = data['structure'])

  restingsets = {}
  for rs_id, data in sstats_dict['resting-sets'].items():
    rs_complexes = [complexes[c_id] for c_id in data['complexes']]
    restingsets[rs_id] = dna.RestingSet(name = data['name'], complexes = rs_complexes)

  reactions = {}
  for rxn_id, data in sstats_dict['reactions'].items():
    reactants = [complexes[c_id] for c_id in data['reactants']]
    products = [complexes[c_id] for c_id in data['products']]
    reactions[rxn_id] = dna.Reaction(name = data['name'],
        reactants = reactants, products = products)

  rs_reactions = {}
  for rsrxn_id, data in sstats_dict['resting-set-reactions'].items():
    reactants = [restingsets[rs_id] for rs_id in data['reactants']]
    products = [restingsets[rs_id] for rs_id in data['products']]
    rs_reactions[rsrxn_id] = dna.RestingSetReaction(name = data['name'],
        reactants = reactants, products = products)

  kparams = sstats_dict['initialization_params']['kinda_params']
  mparams = sstats_dict['initialization_params']['multistrand_params']
  nparams = sstats_dict['initialization_params']['nupack_params']
  pparams = sstats_dict['initialization_params']['peppercorn_params']

  sstats = System(
    complexes = list(complexes.values()),
    restingsets = list(restingsets.values()),
    detailed_reactions = list(reactions.values()),
    condensed_reactions = list(rs_reactions.values()),
    enumeration = False,
    kinda_params = kparams,
    multistrand_params = mparams,
    nupack_params = nparams,
    peppercorn_params = pparams)

  for rs_id, data in sstats_dict['resting-set-stats'].items():
    stats = sstats.get_stats(restingsets[rs_id])
    if stats == None:
      print("Warning: Could not match up stored statistics for {0}.".format(restingsets[rs_id]))
      continue
    nupackjob = stats.get_nupackjob()
    threshold = 0
    for key, val in data.items():
      if key == 'similarity_threshold':
        threshold = val
      elif key == 'c_max':
        stats.c_max = val
      else:
        c = complexes[key]
        c_data = np.array(val['similarity_data'])
        nupackjob.set_complex_prob_data(c.name, c_data)
        if nupackjob.total_sims == 0:
          nupackjob.total_sims = len(c_data)
        else:
          assert nupackjob.total_sims == len(c_data)
    assert threshold > 0
    stats.set_similarity_threshold(threshold)

  for rsrxn_id, data in sstats_dict['resting-set-reaction-stats'].items():
    stats = sstats.get_stats(rs_reactions[rsrxn_id])
    if stats == None:
      print("Warning: Could not match up stored statistics for {0} with a resting-set reaction in the new KinDA object.".format(rs_reactions[rsrxn_id]))
      continue
    multijob = stats.get_multistrandjob()
    stats.multijob_tag = data['tag']

    num_sims = len(data['simulation_data']['tags'])
    sim_data = {key:np.array(d) for key,d in data['simulation_data'].items()}
    multijob.set_simulation_data(sim_data)
    multijob.set_invalid_simulation_data(data['invalid_simulation_data'])
    multijob.total_sims = num_sims

  return sstats

def _import_data_convert_version(sstats_dict, version):
  version_parts = [int(v) for v in version[1:].split('.')]
  if len(version_parts) == 2:
    major, minor = version_parts
    subminor = -1
  elif len(version_parts) == 3:
    major, minor, subminor = version_parts
  else:
    print("KinDA: ERROR: Invalid version number {}. Conversion failed. Simulations and statistical calculations may fail.".format(version))
    return sstats_dict

  if major == 0 and minor == 1 and subminor <= 5:
    print("KinDA: ERROR: Invalid version number {}. Conversion failed. Simulations and statistical calculations may fail.".format(version))
    return sstats_dict
  if major == 0 and minor == 1 and subminor <= 7:
    # add 'valid' entry to all simulation data
    for data in sstats_dict['resting-set-reaction-stats'].values():
      tags = data['simulation_data']['tags']
      num_sims = len(tags)
      MS_TIMEOUT, MS_ERROR = -1, -3
      data['simulation_data']['valid'] = np.array([
        t!=MS_TIMEOUT and t!=MS_ERROR for t in tags])
  if major == 0 and minor == 1 and subminor <= 10:
    # add 'invalid_simulation_data' dict ms_results
    for data in sstats_dict['resting-set-reaction-stats'].values():
      invalid_idxs = [i for i in range(len(data['simulation_data']['valid'])) if data['simulation_data']['valid'][i]==0]
      data['invalid_simulation_data'] = [{'simulation_index': i} for i in invalid_idxs]
  if major == 0 and minor == 1 and subminor <= 12:
    # change macrostate mode name from 'disassoc' to 'ordered-complex'
    if sstats_dict['initialization_params']['kinda_params']['start_macrostate_mode'] == 'disassoc':
      sstats_dict['initialization_params']['kinda_params']['start_macrostate_mode'] = 'ordered-complex'
    if sstats_dict['initialization_params']['kinda_params']['stop_macrostate_mode'] == 'disassoc':
      sstats_dict['initialization_params']['kinda_params']['stop_macrostate_mode'] = 'ordered-complex'

  # sstats_dict should be in current format now
  sstats_dict['version'] = __version__

  return sstats_dict


# Convenience function to create System object from a given PIL file
# Currently only accepts old-style PIL notation (no kernel notation)
def from_pil(path, enumeration = True, **kwargs):
  domains, strands, complexes = io_PIL.from_PIL(path)
  return System(complexes, enumeration = enumeration, **kwargs)


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
            domain_list = [domains[x] for x in line[2]]

            d = dna.Domain(name = name, subdomains = domain_list)
            domains[d.name] = d
            domains[d.complement.name] = d.complement

            # if it is a strand...
            s = dna.Strand(name = name, domains = domain_list)
            strands[s.name] = s
            strands[s.complement.name] = s.complement

            get_strand[tuple(domain_list)] = s


        elif line[0] == 'strand-complex':
            strand_list = [strands[x] for x in line[2]]
            structure = line[3].replace(' ','')
            cplx = dna.Complex(name = name, strands = strand_list, structure = structure )
            complexes[cplx.name] = cplx

        elif line[0] == 'kernel-complex':
            sequence, structure = resolve_loops(line[2])

            # Replace names with domain objects.
            try :
                sequence = [domains[d] for d in sequence]
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
            cplxs = [complexes[c] for c in line[2]]
            resting[name] = dna.RestingSet(name = name, complexes = cplxs)

        elif line[0] == 'reaction':
            rtype = line[1][0][0] if line[1] != [] and line[1][0] != [] else None

            assert rtype is not None

            if rtype == 'condensed' :
                reactants = [resting[c] for c in line[2]]
                products  = [resting[c] for c in line[3]]
                con_reactions.append(
                        dna.RestingSetReaction(reactants = reactants, products = products))
            else :
                reactants = [complexes[c] for c in line[2]]
                products  = [complexes[c] for c in line[3]]
                det_reactions.append(
                        dna.Reaction(reactants = reactants, products = products))

        else :
            print('# Ignoring keyword: {}'.format(line[0]))

    # Make sure the reverse reaction between every pair of reactants is
    # included. These unproductive reactions will be important stop states for
    # Mulstistrand simulations.
    reactant_pairs = it.product(list(resting.values()), list(resting.values()))
    for reactants in reactant_pairs:
        con_reactions.append(dna.RestingSetReaction(reactants = reactants,
                                                    products  = reactants))

    return list(complexes.values()), det_reactions, list(resting.values()), con_reactions

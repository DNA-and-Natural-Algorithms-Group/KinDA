#!/usr/bin/env python3

"""
--------------------------
KinDA executable workflow:
--------------------------

cat input.pil | peppercorn -c -d [--options] | KinDA [--options]

    Selected options:
        --backup kotani.db (exports results to kotani.db)
        --restore kotani.db (loads system from kotani.db)
        --database kotani.db (combines --backup and --restore)

    Potential improvements:
        *) it is not possible to load and extend a database, or to read out
        just the relevant parts. (e.g. reporter systems might be redundant
        between experiments, so why calculate them every time new?)
        Unfortunately, I think one needs to store canonical complex forms for
        that (possible with dsdobjects)

        *) the existing import/export structure is a little more rigid than
        necessary. A clearer separation of system and session-parameters would
        be nice (e.g. multiprocessing)
"""

import os
import sys
import argparse
from typing import List

import peppercornenumerator
from peppercornenumerator.input import read_pil as pepper_read_pil
from peppercornenumerator.output import write_pil as pepper_write_pil

import kinda
from kinda.objects.io_KinDA import read_pil, write_pil, import_data, export_data
from kinda.statistics.stats import RestingSetRxnStats, RestingSetStats
from kinda.simulation.nupackjob import NupackSampleJob
from kinda.simulation.multistrandjob import MultistrandJob


def init_parameter_dicts(args):
    """Initialize kinda's parameter format.
    
    When using the kinda library, parameter defaults are set during
    installation in kinda.options. Here, we want to keep an overview
    and control parameters at runtime.

    Args:
        args (argparse.ArgumentParser()): Commandline parameters.

    Not all parameters can set via commandline arguments, time will tell...
    """
    kparams = {
        # System Parameters
        'start_macrostate_mode': args.macrostate_mode if args.macrostate_mode else \
                                 args.start_macrostate_mode,
        'stop_macrostate_mode': args.macrostate_mode if args.macrostate_mode else \
                                args.stop_macrostate_mode,
        'multistrand_similarity_threshold': args.multistrand_similarity_threshold,
        'nupack_similarity_threshold': args.nupack_similarity_threshold,

        # Session Parameters
        'multistrand_multiprocessing': not args.no_multiprocessing,
        'nupack_multiprocessing': not args.no_multiprocessing,
        'max_concentration': args.max_concentration,
    }
    
    mparams = {
        # Session Parameter
        'verbosity': 0,
        # System Parameters
        'dangles': args.dangles.title(), # first letter uppercase
        'substrate_type': args.material,
        'temperature': args.temperature,
        'simulation_time': args.multistrand_timeout,
        'gt_enable': True, # Always True, to be consistent with NUPACK.
        'sodium': args.sodium,
        'magnesium': args.magnesium,
        'output_interval': -1,
        'parameter_type': 'Nupack',
        'rate_model': 'DNA23Metropolis',
        'join_concentration': 1e-15  # collition rate of complexes. ([] * k_bi)
        # join_concentration must be low enough to ensure that no bind
        # reactions occur after the first step of a simulation.
    }

    nparams = {
        # System Parameters
        'ensemble': args.dangles,
        'material': args.material.lower(),
        'celsius': args.temperature,
        'sodium': args.sodium,
        'magnesium': args.magnesium
    }

    # More Session parameters
    rate_params = {
        'relative_error':  args.error_goal if args.error_goal is not None else \
                           args.rate_error_goal,
        'min_batch_size':  args.rate_batch_size[0],
        'init_batch_size': args.rate_batch_size[1],
        'max_batch_size':  args.rate_batch_size[2],
        'max_sims':        args.max_sims if args.max_sims is not None else \
                           args.rate_max_sims
    }
    rate_params['min_batch_size']  = min([rate_params['min_batch_size'], rate_params['max_sims']])
    rate_params['init_batch_size'] = min([rate_params['init_batch_size'], rate_params['max_sims']])

    prob_params = {
        'relative_error':  args.error_goal if args.error_goal is not None else args.prob_error_goal,
        'min_batch_size':  args.prob_batch_size[0],
        'init_batch_size': args.prob_batch_size[1],
        'max_batch_size':  args.prob_batch_size[2],
        'max_sims':        args.max_sims if args.max_sims is not None else args.prob_max_sims
    }
    prob_params['min_batch_size']  = min([prob_params['min_batch_size'], prob_params['max_sims']])
    prob_params['init_batch_size'] = min([prob_params['init_batch_size'], prob_params['max_sims']])

    return kparams, mparams, nparams, rate_params, prob_params
 

def peppercorn(pilstring, is_file = False, params=None):
    cxs, rxns, strands = pepper_read_pil(pilstring, is_file, composite = True)
    enum = peppercornenumerator.Enumerator(cxs.values(), rxns)
    enum.enumerate()
    return pepper_write_pil(enum, fh = None,
            detailed = True, condensed = True, composite = strands)


def calculate_all_complex_probabilities(KindaSystem, spurious, nsth, multip = True,
        backup = None, verbose = 0, use_pickle = True, **kwargs):
    """ TODO

    Args:
        KindaSystem():
        verbose (int):
        kwargs (dict): Arguments that are passed on to get_conformation_probs() of the
            kinda.statistics.stats.RestingSetRxnStats() object.

    Double check: Probability of the spurious conformation None.
    """
    ## Get resting sets from database, but don't do nothing yet
    restingsets = KindaSystem.get_restingsets(spurious = spurious)

    ## Analyze each resting set
    for e, rms in enumerate(restingsets, 1):
        if verbose:
            print("\n# Analyzing resting set {}/{}: {} with similarity threshold = {}".format(
                e, len(restingsets), rms, nsth))
    
        # Get present stats object
        rms_stats = KindaSystem.get_stats(rms)
        rms_stats.sampler._multiprocessing = multip
        
        if nsth != rms_stats.get_similarity_threshold():
            rms_stats.set_similarity_threshold(nsth)

        # Query/calculate probabilities for all conformations 
        # (including the spurious conformation, denoted as None)
        num = rms_stats.get_num_sims()
        rms_stats.get_conformation_probs(verbose = verbose, spurious = False, **kwargs)

        if num != rms_stats.get_num_sims() and backup:
            export_data(KindaSystem, backup, use_pickle)

        if verbose >= 2:
            tot = 0
            for c in rms.complexes + [None]:
                name = None if c is None else c.name
                tot += rms_stats.get_conformation_prob(name, max_sims=0)
            print("# Sum of probabilities: {}".format(tot))


def calculate_all_reaction_rates(KindaSystem, unproductive, spurious, multip = True,
        backup = None, verbose = 0, use_pickle = True, **kwargs):
    """Calculates reaction rates and error bars.

    There are three types of reactions:
        a) the regular reactions enumerated using peppercorn (and more?)
        b) unproductive reactions

    """
    rxns = KindaSystem.get_reactions(spurious = spurious, unproductive = unproductive)

    ## Simulate each reaction
    for e, rxn in enumerate(rxns, 1):
        if verbose :
            reactants = map(lambda x:x.name, rxn.reactants)
            products  = map(lambda x:x.name, rxn.products)
            print("\n# Analyzing reaction {}/{}: {} -> {}".format(e, len(rxns), 
                ' + '.join(reactants), ' + '.join(products)))
    
        # Get stats object
        rxn_stats = KindaSystem.get_stats(rxn)
        rxn_stats.multijob.multiprocessing = multip
    
        # Query/calculate k1 and k2 reaction rates to requested precision
        goal = kwargs['relative_error']
        maxs = kwargs['max_sims']
        min_bs = kwargs['min_batch_size']
        ini_bs = kwargs['init_batch_size']
        max_bs = kwargs['max_batch_size']

        num = rxn_stats.get_num_sims()
        v = 2 if verbose == 1 else verbose
        rxn_stats.get_raw_stat('k1', goal, maxs, 
                init_batch_size = ini_bs, 
                min_batch_size = min_bs, 
                max_batch_size = max_bs, 
                verbose = v)
        maxs -= (rxn_stats.get_num_sims() - num)
        assert maxs >= 0
        rxn_stats.get_raw_stat('k2', goal, maxs, 
            init_batch_size = ini_bs, 
            min_batch_size = min_bs, 
            max_batch_size = max_bs, 
            verbose = verbose)

        if num != rxn_stats.get_num_sims() and backup:
            export_data(KindaSystem, backup, use_pickle)


def merge_databases(ref_sys: kinda.System, databases: List[str],
                    use_pickle: bool) -> None:
    """
    Merge multiple database files into your current system.

    NOTE: This imports only simulation statistics for resting sets and their
    condensed reactions from the supplied `databases`. For a description of the
    intended use case, see `case_studies/Fig9_Kotani2017/README.md`.
    """
    for db in databases:
        print('# Importing {}'.format(db))
        new_sys = import_data(db, use_pickle)

        ref_rs = {rs.name: rs for rs in ref_sys._restingsets}
        for rs in new_sys._restingsets:
            assert ref_rs[rs.name] == rs
            ref_stats = ref_sys.get_stats(rs)
            if ref_stats is None:
                raise Exception(f"Cannot find resting set: {rs}")
            new_stats = new_sys.get_stats(rs)
            assert isinstance(ref_stats, RestingSetStats)
            assert isinstance(new_stats, RestingSetStats)
            nupackjob = ref_stats.get_nupackjob()
            assert isinstance(nupackjob, NupackSampleJob)

            # Now transfer the data
            new_data = []
            for cx in rs.complexes:
                ref_data = list(nupackjob.get_complex_prob_data(cx.name))
                new_data = list(new_stats.get_conformation_prob_data(cx.name))
                nupackjob.set_complex_prob_data(cx.name, ref_data + new_data)
            nupackjob.total_sims += len(new_data)
            nupackjob.recompute_complex_counts()

        ref_rxns = {repr(rxn): rxn for rxn in ref_sys._condensed_reactions}
        seen_reactants = set()
        for rxn in new_sys._condensed_reactions:
            assert ref_rxns[repr(rxn)] == rxn
            ref_stats = ref_sys.get_stats(rxn)
            if ref_stats == None:
                raise Exception(f"Cannot find resting set reaction: {rxn}")
            new_stats = new_sys.get_stats(rxn)
            assert isinstance(ref_stats, RestingSetRxnStats)
            assert isinstance(new_stats, RestingSetRxnStats)
            assert ref_stats.multijob_tag == new_stats.multijob_tag
            multijob = ref_stats.get_multistrandjob()
            assert isinstance(multijob, MultistrandJob)

            new_reactants = tuple(sorted(rxn.reactants))
            if new_reactants not in seen_reactants:
                seen_reactants.add(new_reactants)

                # Now transfer the data
                multijob.add_simulation_data(new_stats.get_simulation_data())
                multijob.set_invalid_simulation_data(
                    multijob.get_invalid_simulation_data()
                    + new_stats.get_invalid_simulation_data())


def main(args):
    """Calculate DSD system parameters using Mulstistrand and NUPACK.

    KinDA.System is essentially a lazy database. It stores Nupack/Multistrand
    derived parameters from every calculation so far, and updates them whenever
    new values or smaller error bars are requested.

    We distinguish "system parameters", which have to be consistent when
    restoring data from a previously created database file (e.g. temperature,
    material, etc.) and "session parameters", which can be changed freely (e.g.
    error goals, parallelization, etc.).

    # There are multiple points of entry:
    #   1) Use a *.pil file as produced by peppercorn v0.6 or greater.
    #   2) Use a *.pil file and enumerate using KinDA with *installed* default parameters
    #   3) Use --restore file.db to restore a System (ignores STDIN)
    """

    # Input/parameter processing ...
    systeminput = args.input_filename

    kparams, mparams, nparams, rparams, pparams = init_parameter_dicts(args)

    # Import/Export
    if args.database:
        if (args.restore_json or args.backup_json or args.restore or args.backup):
            raise SystemExit('You may not specify any other I/O parameters when using --database.')
        args.restore = args.backup = args.database

    import_pickle = False if args.restore_json else True
    export_pickle = False if args.backup_json else True
    args.backup = args.backup_json if args.backup_json else args.backup
    args.restore = args.restore_json if args.restore_json else args.restore

    # Session parameters
    unproductive = args.unproductive_reactions
    if unproductive is True: 
        unproductive = None

    spurious = args.spurious_reactions
    if spurious is True: 
        spurious = None

    if not args.force and args.backup and os.path.exists(args.backup):
        if args.restore != args.backup:
            raise SystemExit('Use --force to overwrite existing backup file.')

    # System initialization
    if args.restore:
        # Restore a backup file from a previous analysis. This mode ignores any
        # other form of input. Sanity checks ensure that the system parameters
        # match those of the restored database.
        if systeminput :
            raise SystemExit('--restore mode cannot handle *.PIL file input.')
        if not os.path.exists(args.restore):
            raise SystemExit('Cannot find file to restore database: {}.'.format(args.restore))

        if args.verbose:
            print('# Importing options and parameters from {}.'.format(args.restore))

        KindaSystem = import_data(args.restore, import_pickle)

        session_params = set(['nupack_multiprocessing', 'multistrand_multiprocessing',
                'nupack_similarity_threshold', 'max_concentration'])

        for k,v in KindaSystem.initialization_params['kinda_params'].items():
            if k not in kparams:
                print("# WARNING: Using undefined library parameter: {} = {} ".format(k, v))
            elif k in session_params:
                pass
            elif v != kparams[k]:
                print("# WARNING: recovered system uses",
                        "different KinDA parameter: {} = {}".format(k, v))
        for k,v in KindaSystem.initialization_params['multistrand_params'].items():
            if k not in mparams:
                print("# WARNING: Using undefined Multistrand parameter: {} = {} ".format(k, v))
            elif v != mparams[k]:
                print("# WARNING: recovered system uses",
                        "different Multistrand parameter: {} = {}".format(k, v))
        for k,v in KindaSystem.initialization_params['nupack_params'].items():
            if k not in nparams:
                print("# WARNING: Using undefined NUPACK parameter: {} = {} ".format(k, v))
            elif v != nparams[k]:
                print("# WARNING: recovered system uses",
                        "different NUPACK parameter: {} = {}".format(k, v))

    else :
        if not systeminput :
            systeminput = ''
            for l in sys.stdin:
                systeminput += l

        cxs, det, rss, con = read_pil(systeminput, args.input_filename is not None)
        if cxs == []:
            raise SystemExit('ERROR: No complexes found.')
        elif not (det and rss and con): 
            if args.verbose:
                print('# Enumerating using peppercornenumerator-{} with default parameters.'.format(
                    peppercornenumerator.__version__))
            pilstring = peppercorn(systeminput, args.input_filename is not None)
            cxs, det, rss, con = read_pil(pilstring, is_file=False)

        if args.verbose:
            print('# Initializing KinDA system.')
        KindaSystem = kinda.System(cxs, rss, det, con, enumeration = False,
                kinda_params = kparams, multistrand_params = mparams, nupack_params = nparams)

    if args.merge:
        print('# Merging results from given database files:')
        merge_databases(KindaSystem, args.merge, import_pickle)

    # Now that we have the Kinda System setup, we can do two things:
    #   1) calculate probabilities of being in a particular resting complex using NUPACK
    #   2) calculate reaction rates using Multistrand

    # Update complex concentration parameters...
    c_max = dict(term.split('=') for term in args.c_max) if args.c_max else dict()
    if args.verbose:
        print("\n# Setting maximum concentrations")
    restingsets = KindaSystem.get_restingsets(spurious = spurious)
    for x in c_max.keys():
        if x not in map(lambda r: r.name, restingsets):
            print("# Ignoring c-max specification for unknown resting-set: {}".format(x))

    for rms in restingsets:
        rms_stats = KindaSystem.get_stats(rms)
        try:
            rms_stats.c_max = float(c_max[rms.name])
        except KeyError:
            rms_stats.c_max = args.max_concentration
        if args.verbose:
            print("# {} = {} nM".format(rms, 1e9 * rms_stats.c_max))

    # let's do 1)
    calculate_all_complex_probabilities(
        KindaSystem, spurious, args.nupack_similarity_threshold,
        not args.no_multiprocessing, args.backup, args.verbose, export_pickle,
        **pparams)

    # let's do 2)
    calculate_all_reaction_rates(
        KindaSystem, unproductive, spurious, not args.no_multiprocessing,
        args.backup, args.verbose, export_pickle, **rparams)
   
    if args.backup and (not os.path.exists(args.backup) or args.merge):
        export_data(KindaSystem, args.backup, export_pickle)
        if args.verbose:
            print("\n# Results stored in {} using {} format.".format(
                args.backup, 'PICKLE' if export_pickle else 'JSON'))

    ######################
    # Print the results: #
    ######################

    if args.output is None:
        write_pil(KindaSystem, sys.stdout, spurious=spurious, unproductive=unproductive)
    else:
        with open(args.output, 'w') as pil:
            write_pil(KindaSystem, pil, spurious=spurious, unproductive=unproductive)
        print("\n# Results wrote to {} using *.pil format.".format(args.output))


def add_kinda_args(parser):
    interface = parser.add_argument_group('KinDA I/O parameters')
    session   = parser.add_argument_group('KinDA session parameters')
    system    = parser.add_argument_group('KinDA system parameters')

    parser.add_argument('--version', action='version', version='%(prog)s ' + kinda.__version__)
    parser.add_argument( '-v', '--verbose', action='count', default=0,
        help="Print more output (-vv for extra debugging information)")

    # I/O parameters
    interface.add_argument('input_filename', default=None, nargs='?', metavar='<str>',
            help="Path to the input file.")

    interface.add_argument('-o', '--output', default=None, metavar='<str>',
            help="Specify output file. Prints to STDOUT if not provided.")

    interface.add_argument('-b', '--backup', default=None, metavar='<str>',
        help="""Store system in the given pickle file.""")

    interface.add_argument('--backup-json', default=None, metavar='<str>',
        help="""Store system in the given JSON file.""")

    interface.add_argument('-r', '--restore', default=None, metavar='<str>',
        help="""Restore system from given pickle file.""")

    interface.add_argument('--restore-json', default=None, metavar='<str>',
        help="""Restore system from given JSON file.""")

    interface.add_argument('-d', '--database', default=None, metavar='<str>',
        help="""Restore from and backup to given database file (uses PICKLE).""")
        # NOTE: It would be nice to store/load/update *any* given database
        # file, but that requires dsdobjects.

    interface.add_argument('--force', action='store_true',
        help="""Overwrite existing files.""")

    interface.add_argument('--merge', default=None, nargs='+', metavar='<str>',
        help="""Merge a list of database files into your primary analysis
        setup.  Be careful, there (intentionally) no sanity checks for system
        parameters of merged databases. """)

    # Session parameters
    session.add_argument('--unproductive-reactions', action='store_true',
            help="""Toggle to include unproductive reactions for analysis.""")

    session.add_argument('--spurious-reactions', action='store_true',
            help="""Toggle to include spurious reactions and complexes for analysis.""")

    session.add_argument('-e', '--error-goal', type=float, default = None,
            metavar='<float>', 
            help="""If provided, overwrites both --rate-error-goal and
            --prob-error-goal.""")

    session.add_argument('--rate-error-goal', type=float, default = 0.3, 
            metavar='<float>',
            help="""Relative error goal for reaction rates.""")

    session.add_argument('--prob-error-goal', type=float, default = 0.3, 
            metavar='<float>',
            help="""Relative error goal for complex probabilities.""")

    session.add_argument('-s', '--max-sims', type=int, default = None,
            metavar='<int>', 
            help="""If provided, overwrites both --rate-max-sims and
            --prob-max-sims. --max-sims 0 will turn off all simulations.""")

    session.add_argument('--rate-max-sims', type=int, default = 1000, metavar='<int>',
            help="""Maximum number of Multistrand simulations for this Session.""")

    session.add_argument('--prob-max-sims', type=int, default = 1000, metavar='<int>',
            help="""Maximum number of NUPACK simulations for this Session.""")

    session.add_argument('--rate-batch-size', type=int, nargs=3, 
            default = [50, 100, 1000], metavar='<int>',
            help="""Adjust batch size for Multistrand simulations. 

                - First argument: minimum batch size 
                - Second argument: initial batch size 
                - Third argument: maximum batch size 

            When using multiple cores, the batch size is divided by the number
            of cores.  A smaller batch size can lead to slower run time when
            using multiple cores, as all parallel jobs have to finish in order
            to determine whether a new batch gets submitted. A large batch size
            can lead to excess work in order to reach your error goals.""")

    session.add_argument('--prob-batch-size', type=int, nargs=3,
            default = [100, 100, 1000], metavar='<int>',
            help="""Adjust batch size for NUPACK stochastic backtracking. 
            Equivalent behavior to --rate-batch-size. """)

    session.add_argument('--max-concentration', type=float, default = 1e-7, 
            metavar='<float>',
            help="""Maximum concentration of any resting complex in the system [M].""")

    session.add_argument("--c-max", nargs='+', metavar='<str>=<flt>',
            help="""Vector of maximum restingset concentrations.  E.g. \"--c-max
            S=100e-9 C=10e-9\" overwrites the maximum concentration value for
            species "S" and "C". The default value for all species is set with
            --max-concentration. [M]""")

    session.add_argument('--multistrand-timeout', type=float, default = 1.0, 
            metavar='<float>',
            help="""Maximum Multistrand simulation time [seconds].""")

    session.add_argument('--no-multiprocessing', action="store_true",
            help="""Switch off multiprocessing for Multistrand and NUPACK.""")

    session.add_argument('--nupack-similarity-threshold', type=float, 
            default = 0.51, metavar='<float>',
            help="""Calculate complex probabilities (p-approximation) using this 
            similarity threshold. Setting a higher threshold can lead to unexpected 
            results for macrostates that contain more than one complex, as
            (equally likely) transitions between the complexes along branch
            migration domains are not part of this specification.  """)

    # System parameters
    system.add_argument('--macrostate-mode', action="store",
            choices=('ordered-complex', 'count-by-complex', 'count-by-domain'),
            default = None, 
            help="""Use corresponding Multistrand macrostate model. 

                (1) ordered-complex: A complex with strands in specific order.

                (2) count-by-complex: Use p-approximation for complexes. That is,
                  the number of correctly paired and unpaired bases(!) divided by
                  total number of bases.

                (3) count-by-domain: Use p-approximation for each domain. That
                  is, the number of correctly paired and unpaired bases divided
                  by total number of bases in that domain. 
                
            Use --multistrand-similarity-threshold to select if a given structure
            is within the count-by-domain or count-by-complex macrostate.  If
            provided, overwrites both --start-macrostate-mode and
            --stop-macrostate-mode. Defaults are set there. """)

    system.add_argument('--start-macrostate-mode', action="store",
            choices=('ordered-complex', 'count-by-complex', 'count-by-domain'),
            default = 'ordered-complex', 
            help="""Use corresponding Multistrand macrostate definition for
            sampling of initial conformation.""")

    system.add_argument('--stop-macrostate-mode', action="store", 
            choices=('ordered-complex', 'count-by-complex', 'count-by-domain'),
            default = 'ordered-complex', 
            help="""Stop a Multistrand simulation whenever a trajectory reaches
            this macrostate definition. """)

    system.add_argument('--multistrand-similarity-threshold', type=float, 
            default = 0.51, metavar='<float>',
            help="""Similarity threshold to assign complexes to count-by-domain or 
            count-by-complex macrostate-mode.""")

    system.add_argument('-T','--temperature', type=float, default=25, metavar='<float>',
            help="""Simulation temperature in Celsius.""")

    system.add_argument('--sodium', type=float, default=1.0, metavar='<float>',
            help="""Sodium concentration [M].""")

    system.add_argument('--magnesium', type=float, default=0.0, metavar='<float>',
            help="""Magnesium concentration [M].""")

    # Potentially supported parameters.
    system.add_argument('--dangles', default='some',
            choices=('some', 'all', 'none'),
            help="""NUPACK / Multistrand dangle parameter.""")

    system.add_argument('--material', default='DNA',
            choices=('DNA'),  # RNA not supported
            help=argparse.SUPPRESS)


def cli_main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_kinda_args(parser)
    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    cli_main()

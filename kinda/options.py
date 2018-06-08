kinda_params = {
  'stop_macrostate_mode': 'disassoc', # may be 'disassoc', 'count-by-complex', 'count-by-domain'
  'start_macrostate_mode': 'disassoc', # may be 'disassoc', 'count-by-complex', 'count-by-domain'
  'multistrand_similarity_threshold': 0.51,
  'nupack_similarity_threshold': 0.51,
  'multistrand_multiprocessing': True,
  'nupack_multiprocessing': True,
  'enable_unimolecular_reactions': False,
  'unimolecular_k1_scale': 1000,  # any value >= 1000 should be sufficient
  'max_concentration': 1e-7  # Provides a default max concentration for each resting set, used for system-level scores
}

# These defaults are given directly to Multistrand, unless overridden
# while initializing a System object
multistrand_params = {
  'verbosity': 0,
  'dangles': 'Some',
  'gt_enable': True,
  'simulation_time': 1.00, 
  'output_interval': 0,
  'parameter_type': 'Nupack',
  'substrate_type': 'DNA',
  'rate_method': 'Metropolis',
  'temperature': 25,
  'sodium': 1.0,
  'magnesium': 0.0,
  'unimolecular_scaling': 2.41686715e6,
  'bimolecular_scaling': 8.01171383e5,
  'join_concentration': 1e-15 # this must be low enough to ensure no bind reactions occur after the first step of a simulation
}

# These defaults are given directly to the Nupack Python interface, unless
# overridden while initializing a System object
nupack_params = {
  'dangles': 'some',
  'material': 'dna',
  'T': 25,
  'sodium': 1.0,
  'magnesium': 0.0
}

# These defaults are given directly to Peppercorn, unless overridden
# while initializing a System object
peppercorn_params = {
  'max_complex_size': 6,
  'max_reaction_count': 1000,
  'max_complex_count': 200,
  'release_cutoff_1_1': 8,
  'release_cutoff_1_N': 8,
  'remote_migration': True,
  'max_helix_migration': True
}

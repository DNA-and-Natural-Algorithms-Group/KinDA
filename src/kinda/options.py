kinda_params = {
  # may be 'ordered-complex', 'count-by-complex', 'count-by-domain'
  'stop_macrostate_mode': 'ordered-complex',
  'start_macrostate_mode': 'ordered-complex',

  'multistrand_similarity_threshold': 0.51,
  'nupack_similarity_threshold': 0.51,
  'multistrand_multiprocessing': True,
  'nupack_multiprocessing': True,
  'enable_unimolecular_reactions': False,
  # any value >= 1000 should be sufficient
  'unimolecular_k1_scale': 1000,
  # Provides a default max concentration for each resting set, used for
  # system-level scores
  'max_concentration': 1e-7
}

# These defaults are given directly to Multistrand, unless overridden while
# initializing a System object
multistrand_params = {
  'verbosity': 0,

  'dangles': 'Some',
  'gt_enable': True,

  'simulation_time': 1.00,
  'output_interval': -1,

  'parameter_type': 'Nupack',
  'substrate_type': 'DNA',

  'temperature': 25,
  'sodium': 1.0,
  'magnesium': 0.0,
  # this must be low enough to ensure no bind reactions occur after the first
  # step of a simulation
  'join_concentration': 1e-15,

  # Choose a kinetic model from the presets defined in
  # `multistrand.options.Options`:
  #     "JSDefault", "JSMetropolis25", "JSKawasaki25", "JSKawasaki37",
  #     "JSMetropolis37", "DNA23Metropolis", "DNA23Arrhenius", "DNA29Arrhenius"
  'rate_model': "DNA23Metropolis",

  # Alternatively, use the following options to specify a custom kinetic model.
  # 'rate_method': 'Metropolis',
  # 'unimolecular_scaling': 2.41686715e6,
  # 'bimolecular_scaling': 8.01171383e5,
}

# These defaults are given directly to the Nupack Python interface, unless
# overridden while initializing a System object
nupack_params = {
  'ensemble': 'some',
  'material': 'dna',
  'celsius': 25,
  'sodium': 1.0,
  'magnesium': 0.0
}

# These defaults are given directly to Peppercorn, unless overridden while
# initializing a System object
peppercorn_params = {
  'max_complex_size': 6,
  'max_reaction_count': 1000,
  'max_complex_count': 200,
  'release_cutoff_1_1': 8,
  'release_cutoff_1_N': 8,
  'reject_remote': False,
  'max_helix': True
}

flags = {
  'loose_complexes': False
}

general_params = {
  'loose_complex_similarity': 0.7,
  'material': 'DNA'
}

multistrand_params = {
  'multiprocessing': True,  # set to True to use multiple cores when running simulations
  'options': {  # These options are given directly to Multistrand
    'dangles': 'Some',
    'simulation_time': 1.00, 
    'output_interval': 100,
    'parameter_type': 'Nupack',
    'substrate_type': 'DNA',
    'rate_method': 'Metropolis',
    'temperature': 25,
    'unimolecular_scaling': 5.0e6,
    'bimolecular_scaling': 1.4e6,
    'join_concentration': 1e-15
  }
}

nupack_params = {
  'dangles': 'some',
  'material': 'dna',
  'temp': 25
}

peppercorn_params = {
  '--release-cutoff-1-1': 6,
  '--release-cutoff-1-n': 6
}

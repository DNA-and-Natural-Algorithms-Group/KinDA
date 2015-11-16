flags = {
  'no-cleanup': False,
  'loose_complexes': False
}

general_params = {
  'loose_complex_similarity': 0.7
}

multistrand_params = {
  'dangles': 'Some',
  'sim_time': 1.00, 
  'output_interval': 100,
  'param_type': 'Nupack',
  'substrate_type': 'DNA',
  'rate_method': 'Metropolis',
  'temp': 25,
  'bimolecular_scaling': 0.5e5,
  'multithreading': True,
  'join_concentration': 1e-15
}

nupack_params = {
  'dangles': 'some',
  'material': 'dna',
  'temp': 25
}

peppercorn_params = {
  '--release-cutoff-1-1': 5,
  '--release-cutoff-1-n': 6
}

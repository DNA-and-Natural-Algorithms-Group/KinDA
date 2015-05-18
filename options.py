flags = {
  'no-cleanup': False,
  'loose_complexes': False
}

general_params = {
  'loose_complex_similarity': 0.66
}

multistrand_params = {
  'dangles': 'Some',
  'sim_time': 1e-6,
  'output_interval': 100,
  'param_type': 'Nupack',
  'substrate_type': 'DNA',
  'rate_method': 'Metropolis',
  'temp': 25,
  'bimolecular_scaling': 0.5e5,
  'multithreading': False,
  'join_concentration': 1e-15 # setting this to 0 makes any non-single-stranded complexes have "inf" energy, which could be problematic...
}

nupack_params = {
  'dangles': 'some',
  'material': 'dna',
  'temp': 25
}


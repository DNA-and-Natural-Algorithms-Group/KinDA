import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import kinda


def plot_reaction_data(rxns_reactants, datasets, dataset_labels, dataset_colors, title, figsize):
  num_datasets = len(datasets)
  num_rxns = len(rxns_reactants)

  k1 = np.empty((num_rxns, num_datasets))
  k1_err = np.empty((num_rxns, num_datasets))
  k2 = np.empty((num_rxns, num_datasets))
  k2_err = np.empty((num_rxns, num_datasets))

  plt.figure(figsize = figsize)
  axes_k1 = plt.subplot(211)
  axes_k2 = plt.subplot(212)
  for i,sstats in enumerate(datasets):
    reactants = [[sstats.get_restingset(name = n) for n in reactant_names] for reactant_names in rxns_reactants]
    
    rxns = [sstats.get_reaction(reactants = r, unproductive=False, spurious=False) for r in reactants]
    rxns_stats = [sstats.get_stats(rxn) for rxn in rxns]
  
    k1[:, i] = [rxn_stats.get_k1(max_sims=0) for rxn_stats in rxns_stats]
    k1_err[:, i] = [rxn_stats.get_k1_error(max_sims=0) for rxn_stats in rxns_stats]
    k2[:, i] = [rxn_stats.get_k2(max_sims=0) for rxn_stats in rxns_stats]
    k2_err[:, i] = [rxn_stats.get_k2_error(max_sims=0) for rxn_stats in rxns_stats]
  
    left = np.arange(num_rxns)*(num_datasets+1) + i
    axes_k1.bar(left = left, height = k1[:,i], width=1, yerr = k1_err[:,i], capsize = 5, color = dataset_colors[i], edgecolor = '0.0')
    axes_k2.bar(left = left, height = k2[:,i], width=1, yerr = k2_err[:,i], capsize = 5, color = dataset_colors[i], edgecolor = '0.0')
  
  xticks = np.arange(num_rxns)*(num_datasets+1) + (num_datasets-1)/2.
  xlabels = [' + '.join(r) for r in rxns_reactants]

  axes_k1.legend(dataset_labels)
  axes_k1.set_ylabel('$k_1, M^{-1}s^{-1}$')
  axes_k1.set_xticks(xticks)
  axes_k1.set_xticklabels(xlabels)
  axes_k1.set_yscale('log')
  axes_k1.set_ylim(936816.69397475501, 14643556.512799658)

  axes_k2.legend(dataset_labels)
  axes_k2.set_ylabel('$k_2, s^{-1}$')
  axes_k2.set_xticks(xticks)
  axes_k2.set_xticklabels(xlabels)
  axes_k2.set_yscale('log')
  axes_k2.set_ylim(30.510499687666091, 21417005.100912001)

  axes_k1.set_title(title)

  plt.savefig('{}.pdf'.format(title))

  return k1, k1_err, k2, k2_err

# Change these filenames to use your own data
sstats_OR_disassoc = kinda.import_data('Groves2016_OR_ordered-complex.kinda')
sstats_OR_cbc = kinda.import_data('Groves2016_OR_count-by-complex.kinda')
sstats_OR_cbd = kinda.import_data('Groves2016_OR_count-by-domain.kinda')
sstats_AND_disassoc = kinda.import_data('Groves2016_AND_ordered-complex.kinda')
sstats_AND_cbc = kinda.import_data('Groves2016_AND_count-by-complex.kinda')
sstats_AND_cbd = kinda.import_data('Groves2016_AND_count-by-domain.kinda')

dataset_labels = ['ordered-complex', 'count-by-complex', 'count-by-domain']
dataset_colors = ['1.0', '0.8', '0.6']
OR_k1, OR_k1_err, OR_k2, OR_k2_err = plot_reaction_data([['InA', 'OR'], ['InB', 'OR']], [sstats_OR_disassoc, sstats_OR_cbc, sstats_OR_cbd], dataset_labels, dataset_colors, 'fig10D_OR', figsize=(3.6,4.8))
AND_k1, AND_k1_err, AND_k2, AND_k2_err = plot_reaction_data([['InB', 'AND'], ['InA', 'AND_InB']], [sstats_AND_disassoc, sstats_AND_cbc, sstats_AND_cbd], dataset_labels, dataset_colors, 'fig10D_AND', figsize=(3.6,4.8))

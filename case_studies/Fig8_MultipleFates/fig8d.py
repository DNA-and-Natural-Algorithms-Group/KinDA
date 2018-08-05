# Make plots for Figure 8D (Multiple fates case study)
import numpy as np

import kinda

# Replace these paths with your own data output from fig8_simulate.py if desired.
DATA_PATH1 = 'fig8_unmod_raw.kinda' # path of unmodified system data
DATA_PATH2 = 'fig8_mod_raw.kinda' # path of modified system data

IMGS_PREFIX = 'fig8D'

## Import collected data
sstats1 = kinda.import_data(DATA_PATH1)
sstats2 = kinda.import_data(DATA_PATH2)

## Get the relevant resting sets
reaction_names = ['Gate + Interloper -> Fate1', 'Gate + Interloper -> Fate2']
reactant_names = ['Gate1', 'Interloper']
product_names = [['Fate1_Cpx1', 'Fate1_Cpx2'], ['Fate2_Cpx1', 'Fate2_Cpx2a']]
gate_complex_names = ['Gate1', 'Gate1_open', 'Gate2', 'Gate2_open']

all_k1 = np.empty(shape = (2, len(reaction_names)))
all_k1_err = np.empty(shape = (2, len(reaction_names)))
all_k2 = np.empty(shape = (2, len(reaction_names)))
all_k2_err = np.empty(shape = (2, len(reaction_names)))
all_gate_prob = np.empty(shape = (2, len(gate_complex_names)+1))
all_gate_prob_err = np.empty(shape = (2, len(gate_complex_names)+1))

## Extract reaction rate data
for j, (rxn_name, prod_names) in enumerate(zip(reaction_names,product_names)):
  for i, sstats in enumerate([sstats1, sstats2]):
    reactants = [sstats.get_restingset(complex_name = n) for n in reactant_names]
    products = [sstats.get_restingset(complex_name = n) for n in prod_names]
    rxn_stats = sstats.get_stats(sstats.get_reaction(reactants = reactants, products = products))
    all_k1[i,j] = rxn_stats.get_k1(max_sims=0)
    all_k1_err[i,j] = rxn_stats.get_k1_error(max_sims=0)
    all_k2[i,j] = rxn_stats.get_k2(max_sims=0)
    all_k2_err[i,j] = rxn_stats.get_k2_error(max_sims=0)

## Extract 'Gate' resting set data
for i, sstats in enumerate([sstats1, sstats2]):
  gate = sstats.get_restingset(complex_name = 'Gate1')
  rs_stats = sstats.get_stats(gate)
  for j, complex_name in enumerate(gate_complex_names):
    all_gate_prob[i,j] = rs_stats.get_conformation_prob(complex_name, max_sims=0)
    all_gate_prob_err[i,j] = rs_stats.get_conformation_prob_error(complex_name, max_sims=0)
  all_gate_prob[i,4] = rs_stats.get_conformation_prob(None, max_sims=0)
  all_gate_prob_err[i,4] = rs_stats.get_conformation_prob_error(None, max_sims=0)

import matplotlib.pyplot as plt
plt.ion()

## Figure 8D-left (k1 rate comparison)
plt.figure(figsize = (5.2,4.8))
plt.bar(left = [-0.5, 2], height = all_k1[:,0], width = 1, yerr = all_k1_err[:,0], capsize = 15, color = (.2,.2,.6))
plt.bar(left = [0.5, 3], height = all_k1[:,1], width = 1, yerr = all_k1_err[:,1], capsize=15, color = (.4,.8,.4))
plt.legend(['$\mathcal{R}_1$', '$\mathcal{R}_2$'])
plt.xticks([0,2.5], ['$s$ = ATATAT', '$r$ = GCGCGC'])
plt.ylabel('Second-order rate estimate $\hat{k}_1, M^{-1}s^{-1}$')
plt.yscale('log')
plt.xlim(-2.5, 5)
plt.savefig('{}_k1.pdf'.format(IMGS_PREFIX))

## Figure 8D-right (conformation probability comparison)
colors = [(.2,.2,.6), (.4,.8,.4), (.9,.9,.9)]
left = [1,2]
height = [all_gate_prob[:,0:2].sum(axis=1), all_gate_prob[:,2:4].sum(axis=1), all_gate_prob[:,4]]
bottom = [[0, 0], all_gate_prob[:,0:2].sum(axis=1), all_gate_prob[:,0:4].sum(axis=1)]
yerr = [all_gate_prob_err[:,0:2].sum(axis=1), all_gate_prob_err[:,2:4].sum(axis=1), all_gate_prob_err[:,4]]
plt.figure(figsize = (5.2,4.8))
for h,b,e,c in zip(height, bottom, yerr, colors):
  print h,b,e,c
  plt.bar(left, height=h, bottom=b, yerr=e, color=c, capsize=15, width=0.5)

plt.legend(['$B_1 + B_2$', '$B_3 + B_4$', 'spur.'])
plt.xticks([1,2], ['$s$=ATATAT', '$r$=GCGCGC'])
plt.ylabel('probability')
plt.xlim(.25, 2.75)
plt.savefig('{}_probs.pdf'.format(IMGS_PREFIX))

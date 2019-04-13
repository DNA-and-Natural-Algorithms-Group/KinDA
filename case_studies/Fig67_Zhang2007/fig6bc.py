# Generate Figure 6B-C

import datetime

import numpy as np

import kinda

# Modify the following line to use your own simulation data.
DATA_PATH = 'fig6_raw.kinda'

# Helper function to get temporary depletion matrix
def get_temp_depletion(sstats, restingset_names, restingset_concs):
  ## Retrieve RestingSet objects
  restingsets = [sstats.get_restingset(name=name) for name in restingset_names]

  ## Check resting set concentrations
  for i,rs in enumerate(restingsets):
    rs_stats = sstats.get_stats(rs)
    rs_stats.c_max = restingset_concs[i]
  
  ## In a numpy array, store temporary depletion levels for each pairwise combination of resting sets.
  temp_depletion = np.empty(shape=(len(restingsets), len(restingsets)+1))
  for i,rs1 in enumerate(restingsets):
    rs_stats = sstats.get_stats(rs1)
    unprod_rxns = [sstats.get_reaction(unproductive=True, reactants=[rs1,rs2]) for rs2 in restingsets]
    rxn_stats = [sstats.get_stats(rxn) for rxn in unprod_rxns]
    temp_deps = [rs_stats.get_temporary_depletion_due_to(s, max_sims=0) for s in rxn_stats]
  
    temp_depletion[i][:-1] = temp_deps
    temp_depletion[i][-1] = rs_stats.get_temporary_depletion(max_sims=0) # total depletion

  return temp_depletion
 
## Import collected data
sstats = kinda.import_data(DATA_PATH)

## Settings
restingset_names = ['Substrate', 'Catalyst', 'Signal', 'Intermediate', 'Fuel', 'Output', 'Waste']
restingset_concs1 = [10e-9, 1e-9, 10e-9, 10e-9, 13e-9, 10e-9, 10e-9]
restingset_concs2 = [2.5e-7]*7

temp_depletion_low = get_temp_depletion(sstats, restingset_names, restingset_concs1)
temp_depletion_high = get_temp_depletion(sstats, restingset_names, restingset_concs2)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

names_caps = [s.upper() for s in restingset_names]
VMIN = 0.0
VMAX = 0.30
cbar_ticks = np.sqrt(np.arange(0,VMAX+.005,.01))
cbar_tick_labels = ['{:.1f}%'.format(v) if v%5==0 else '' for v in np.arange(0,100*VMAX+.5,1)]

## Plot temporary depletion for low concentrations (x = 10 nM)
plt.figure()
plt.imshow(np.sqrt(temp_depletion_low), interpolation='nearest', cmap=matplotlib.cm.gray_r, vmin=np.sqrt(VMIN), vmax=np.sqrt(VMAX))
plt.xticks(np.arange(len(names_caps)+1), names_caps + ['TOTAL'], rotation='vertical')
plt.yticks(np.arange(len(names_caps)), names_caps)
cbar = plt.colorbar()
cbar.set_ticks(cbar_ticks)
cbar.ax.set_yticklabels(cbar_tick_labels)
plt.title('x = 10nM')
plt.savefig('fig6b.pdf')

plt.figure()
plt.imshow(np.sqrt(temp_depletion_high), interpolation='nearest', cmap=matplotlib.cm.gray_r, vmin=np.sqrt(VMIN), vmax=np.sqrt(VMAX))
plt.xticks(np.arange(len(names_caps)+1), names_caps + ['TOTAL'], rotation='vertical')
plt.yticks(np.arange(len(names_caps)), names_caps)
cbar = plt.colorbar()
cbar.set_ticks(cbar_ticks)
cbar.ax.set_yticklabels(cbar_tick_labels)
plt.title('x = 250nM')
plt.savefig('fig6c.pdf')

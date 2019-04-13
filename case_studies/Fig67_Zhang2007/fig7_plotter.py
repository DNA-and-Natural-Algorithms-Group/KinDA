# Plot analyzed information for Figure 7 (Entropy-driven catalyst, Zhang et al, Science 2007)
# -- Usage is parallel to fig7_analyze.py, but automatically reads all temperatures, and requires a directory name.

# Usage:
#   python fig7_plotter.py . demo
#   python fig7_plotter.py publication_data randomR

import sys, os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 
font = {'weight' : 'bold', 'size'   : 18}
matplotlib.rc('font', **font)
matplotlib.rcParams['axes.linewidth'] = 2
matplotlib.rcParams['lines.linewidth'] = 3

if not len(sys.argv) == 3: 
  sys.exit('ERROR: Needs 2 arguments: a directory name, and simulation type.')

DIR = sys.argv[1]
MODE = sys.argv[2]

MODELIST = [ 'demo', 'trial', 'publication', 'randomT', 'randomR' ]
if not MODE in MODELIST:
  sys.exit('Unrecognized mode.')

if not DIR[-1] == '/':
  DIR = DIR + '/'

temps=[]

complex_names=[
  '[Complex(Catalyst)]',
  '[Complex(Substrate)]',
  '[Complex(Fuel)]',
  '[Complex(Intermediate)]',
  '[Complex(Output)]',
  '[Complex(Signal)]',
  '[Complex(Waste)]'
  ]
depletion=np.zeros([100,len(complex_names)])
formation=np.zeros([100,len(complex_names),3,3])
[FORM_PROB,FORM_ERR,NUPACK_P] = [0,1,2]
papproximation=np.zeros([100,len(complex_names),25])
short_complex_names=['INPUT (C)','SUBSTRATE (S)','FUEL (F)','INTERMEDIATE (I)','OUTPUT (OB)','SIGNAL (SB)','WASTE (W)']

reaction_names=[
  '[Complex(Catalyst)] + [Complex(Substrate)] -> [Complex(Intermediate)] + [Complex(Signal)]',
  '[Complex(Intermediate)] + [Complex(Signal)] -> [Complex(Catalyst)] + [Complex(Substrate)]',
  '[Complex(Fuel)] + [Complex(Intermediate)] -> [Complex(Catalyst)] + [Complex(Output)] + [Complex(Waste)]',
  '[Complex(Fuel)] + [Complex(Substrate)] -> [Complex(cpx_F:LB:OB)] + [Complex(cpx_SB)]',
  '[Complex(Fuel)] + [Complex(Substrate)] -> [Complex(cpx_F:SB:LB)] + [Complex(cpx_OB)]'
  ]
rates=np.zeros([100,len(reaction_names),4])  
[RXN_K1,RXN_K1_ERR,RXN_K2,RXN_K2_ERR] = [0,1,2,3]
histdata=np.zeros([100,len(reaction_names),30])
histbins=np.zeros([100,len(reaction_names),31])
short_reaction_names=[
  r'C+S $\rightarrow$ I+SB',
  r'I+SB $\rightarrow$ C+S',
  r'F+I $\rightarrow$ C+OB+W',
  r'F+S $\rightarrow$ LEAK1+SB',
  r'F+S $\rightarrow$ LEAK2+OB'
  ]

# first just the sequences

fig = plt.figure(frameon=False)
ax = fig.add_axes([0.1, .5, 0.9, 0.8])
ax.axis('off')
if MODE == 'publication':
  plt.text(0.1, 0.9, "Zhang et al 2007", color='k', fontsize=30)
  plt.text(0.1, 0.8, "d1 = CTTTCCTACA", fontname='monospace', color='b', fontsize=25)
  plt.text(0.1, 0.7, "d2 = CCTACGTCTCCAACTAACTTACGG", fontname='monospace', color='b', fontsize=25)
  plt.text(0.1, 0.6, "t3 = CCCT", fontname='monospace', color='b', fontsize=25)
  plt.text(0.1, 0.5, "d4 = CATTCAATACCCTACG", fontname='monospace', color='b', fontsize=25)
  plt.text(0.1, 0.4, "t5 = TCTCCA", fontname='monospace', color='b', fontsize=25)
  plt.text(0.1, 0.3, "d6 = CCACATACATCATATT", fontname='monospace', color='b', fontsize=25)
elif MODE == 'randomR':
  plt.text(0.1, 0.9, "Modified long domains", color='k', fontsize=30)
  plt.text(0.1, 0.8, "d1 = AACCTGTCCG", fontname='monospace', color='r', fontsize=25)
  plt.text(0.1, 0.7, "d2 = ATGAATCGTAACGTTTTGCAATGG", fontname='monospace', color='r', fontsize=25)
  plt.text(0.1, 0.6, "t3 = CCCT", fontname='monospace', color='b', fontsize=25)
  plt.text(0.1, 0.5, "d4 = TACGGACCTTTAGCGA", fontname='monospace', color='r', fontsize=25)
  plt.text(0.1, 0.4, "t5 = TCTCCA", fontname='monospace', color='b', fontsize=25)
  plt.text(0.1, 0.3, "d6 = GTTGTATAGGCGCAAT", fontname='monospace', color='r', fontsize=25)
elif MODE == 'randomT':
  plt.text(0.1, 0.9, "Modified short domains", color='k', fontsize=30)
  plt.text(0.1, 0.8, "d1 = CTTTCCTACA", fontname='monospace', color='b', fontsize=25)
  plt.text(0.1, 0.7, "d2 = CCTACGTCTCCAACTAACTTACGG", fontname='monospace', color='b', fontsize=25)
  plt.text(0.1, 0.6, "t3 = CGCG", fontname='monospace', color='m', fontsize=25)
  plt.text(0.1, 0.5, "d4 = CATTCAATACCCTACG", fontname='monospace', color='b', fontsize=25)
  plt.text(0.1, 0.4, "t5 = TATTAA", fontname='monospace', color='m', fontsize=25)
  plt.text(0.1, 0.3, "d6 = CCACATACATCATATT", fontname='monospace', color='b', fontsize=25)
plt.savefig("plot0.pdf".format(MODE), bbox_inches='tight')
plt.close()

# now the actual data plots

for temp in range(0,100):
  ANALYSIS_PATH = '{}figT{:d}_{}_analysis.csv'.format(DIR,temp,MODE)
  if os.path.exists(ANALYSIS_PATH):
    print "Loading data for temperature T = {:d} and {} mode.".format(temp,MODE)
    temps.append(temp)

    num_p=0
    with open(ANALYSIS_PATH, 'r') as f:
      reader = csv.reader(f)
      for row in reader:
        if len(row)>0 and "=" not in row[0] and "RESTING" in row[0]:
          nupack_p = { 'dem':0.90, 'tri':0.90, 'pub':0.70, 'ran':0.70 }[MODE[0:3]] 
          num_p = 1
        if len(row)>0 and "=" in row[0] and "RESTING" in row[0]:
          i = row[0].index('=')
          nupack_p = float(row[0][i+2:len(row[0])])
          num_p = num_p + 1
        if len(row)==8 and not row[0].find('->') == -1 and row[0] in reaction_names:
          rxn_i = reaction_names.index(row[0])
          rates[temp,rxn_i] = [ float(row[1]), float(row[2]), float(row[3]), float(row[4]) ]
        if len(row)==4 and not row[0].find('Complex') == -1 and row[0].find('->') == -1:  
          cmplx_i = complex_names.index(row[0])
          depletion[temp,cmplx_i] = float(row[3]) 
          formation[temp,cmplx_i,num_p-1] = [ float(row[1]), float(row[2]), nupack_p ]
        if len(row)==26 and not row[0].find('Complex') == -1 and row[0].find('->') == -1:  
          cmplx_i = complex_names.index(row[0])
          for v in range(25): # p = 0.51, 0.53, ..., 0.97, 0.99
            papproximation[temp,cmplx_i,v] = float(row[v+1])
        if len(row)==62 and not row[0].find('->') == -1 and row[0] in reaction_names:
          rxn_i = reaction_names.index(row[0])
          for h in range(30):
            histdata[temp,rxn_i,h] = float(row[h+1])
          for b in range(31):
            histbins[temp,rxn_i,b] = float(row[b+31])

print num_p

print

## double-check that some stuff got read correctly, if you want
if False:  
  for temp in temps:
    for i in range(len(reaction_names)):
      print "T = {} : rxn{}_k1 = {} +/- {} /M/s, rxn{}_k2 = {} +/- {} /s".format(temp,i,rates[temp,i,0],rates[temp,i,1],i,rates[temp,i,2],rates[temp,i,3])

  print

  for temp in temps:
    for i in range(len(complex_names)):
      for pi in range(num_p):
        print "T = {}, {} : {}-well-formed = {:4.1f} +/- {:4.1f} %, temporary depletion = {:4.1f} %".format(temp,
          complex_names[i],formation[temp,i,pi,NUPACK_P],100*formation[temp,i,pi,FORM_PROB],100*formation[temp,i,pi,FORM_ERR],100*depletion[temp,i])

colors=['r','b','g','m','c','k','y']

all_k1 = [ rates[temp,i,RXN_K1] for temp in temps for i in range(len(reaction_names)) ]
axisminmax = [ min(temps)-5, max(temps)+5, 0.005*min(all_k1), 2*max(all_k1) ]  # [xmin xmax ymin ymax]
axisminmax = [ 10, 65, 3, 9 * 10**6 ]  # for consistency

plt.figure()
for rxn_i in range(len(reaction_names)):
  t  = [ temp + rxn_i / 10.0 for temp in temps if not np.isnan(rates[temp,rxn_i,RXN_K2]) ]
  k1 = [ rates[temp,rxn_i,RXN_K1] for temp in temps if not np.isnan(rates[temp,rxn_i,RXN_K2]) ]
  k1e = [ rates[temp,rxn_i,RXN_K1_ERR] for temp in temps if not np.isnan(rates[temp,rxn_i,RXN_K2]) ]
  plt.errorbar(t,k1,yerr=k1e,fmt=colors[rxn_i]+'o-')

if False: ## do we want to show the k2 bounds when simulations had no successful trajectories?
  for rxn_i in range(len(reaction_names)):
    t  = [ temp + rxn_i / 10.0 for temp in temps if np.isnan(rates[temp,rxn_i,RXN_K2]) ]
    k1 = [ rates[temp,rxn_i,RXN_K1] for temp in temps if np.isnan(rates[temp,rxn_i,RXN_K2]) ]
    k1e = [ rates[temp,rxn_i,RXN_K1_ERR] for temp in temps if np.isnan(rates[temp,rxn_i,RXN_K2]) ]
    plt.errorbar(t,k1,yerr=k1e,fmt=colors[rxn_i]+'s')

# plt.title(MODE + ':     ' + 'Bimolecular rate constants')  
plt.yscale('log')
plt.axis(axisminmax)
plt.xlabel('temperature, C', weight='bold')
plt.ylabel('reaction k$_1$, M$^{-1}$s$^{-1}$', weight='bold')
# plt.legend(short_reaction_names,loc='lower right')
plt.savefig("plot1.pdf", bbox_inches='tight')
plt.close()

all_k2 = [ rates[temp,i,RXN_K2] for temp in temps for i in range(len(reaction_names)) ]
axisminmax = [ min(temps)-5, max(temps)+5, 0.1*min(all_k2), 2*max(all_k2) ]  # [xmin xmax ymin ymax]
axisminmax = [ 10, 65, 0.05, 2 * 10**4 ]  # for consistency

plt.figure()
for rxn_i in range(len(reaction_names)):
  t  = [ temp + rxn_i / 10.0 for temp in temps if not np.isnan(rates[temp,rxn_i,RXN_K2]) ]
  k2 = [ rates[temp,rxn_i,RXN_K2] for temp in temps if not np.isnan(rates[temp,rxn_i,RXN_K2]) ]
  k2e = [ rates[temp,rxn_i,RXN_K2_ERR] for temp in temps if not np.isnan(rates[temp,rxn_i,RXN_K2]) ]
  plt.errorbar(t,k2,yerr=k2e,fmt=colors[rxn_i]+'o-')

# plt.title('Unimolecular rate constants')  
plt.yscale('log')
plt.axis(axisminmax)
plt.xlabel('temperature, C', weight='bold')
plt.ylabel('reaction k$_2$, s$^{-1}$', weight='bold')
if MODE == "randomT":
  plt.legend(short_reaction_names,loc='lower right')
plt.savefig("plot2.pdf", bbox_inches='tight')
plt.close()

all_dep = [ depletion[temp,i] for temp in temps for i in range(len(complex_names))]
axisminmax = [ min(temps)-5, max(temps)+5, -0.05, 1.05 ]  # [xmin xmax ymin ymax]
axisminmax = [ 10, 65, -0.05, 1.05 ]  # for consistency

plt.figure()
for cmplx_i in range(len(complex_names)):
  t  = [ ti + cmplx_i / 10.0 for ti in temps ]
  dep = [ depletion[temp,cmplx_i] for temp in temps ]
  plt.plot(t,dep,colors[cmplx_i]+'o-')

# plt.title('Temporary depletion for 100 nM concentrations')  
plt.axis(axisminmax)
plt.xlabel('temperature, C', weight='bold')
plt.ylabel('fraction temporary depletion', weight='bold')
if MODE == "randomR":
  plt.legend(short_complex_names,loc='upper right')
plt.savefig("plot3.pdf", bbox_inches='tight')
plt.close()

all_form = [ formation[temp,i,0] for temp in temps for i in range(len(complex_names))]
axisminmax = [ min(temps)-5, max(temps)+5, -0.05, 1.05 ]  # [xmin xmax ymin ymax]
axisminmax = [ 10, 65, -0.05, 1.05 ]  # for consistency


for pi in range(num_p):
  nupack_p = formation[temp,0,pi,NUPACK_P]

  plt.figure()
  for cmplx_i in range(len(complex_names)):
    t  = [ ti + cmplx_i / 10.0 for ti in temps ]
    form  = [ formation[temp,cmplx_i,pi,FORM_PROB] for temp in temps ]
    forme = [ formation[temp,cmplx_i,pi,FORM_ERR] for temp in temps ]
    plt.errorbar(t,form,yerr=forme,fmt=colors[cmplx_i]+'o-')

  # plt.title('Fraction of {}-approximate conformations'.format(nupack_p))  
  plt.axis(axisminmax)
  plt.xlabel('temperature, C', weight='bold')
  plt.ylabel('fraction well-formed', weight='bold')
  # plt.legend(short_complex_names,loc='lower right')
  plt.savefig("plot{}.pdf".format(pi+4), bbox_inches='tight')
  plt.close()

if num_p==1:
  cmd="pdfjam plot0.pdf plot1.pdf plot2.pdf plot3.pdf plot4.pdf --nup 1x5 --outfile figTEMP_"+MODE+".pdf --papersize '{1.4in,4.5in}'"
  os.system(cmd)
  cmd="rm plot0.pdf plot1.pdf plot2.pdf plot3.pdf plot4.pdf"
  os.system(cmd)
elif num_p==3:
  cmd="pdfjam plot1.pdf plot2.pdf plot3.pdf plot4.pdf plot5.pdf plot6.pdf --nup 3x2 --outfile figTEMP_"+MODE+".pdf --papersize '{6in,3in}'"
  os.system(cmd)
  cmd="rm plot1.pdf plot2.pdf plot3.pdf plot4.pdf plot5.pdf plot6.pdf"
  os.system(cmd)
else:
  print('Not sure how to arrange plots.  Did not combine plot*.pdf.')

### let's see how p affects the fraction of p-approximations, just for T=25

plt.figure()
pv = [round(0.02*v + 0.51,2) for v in range(25)]
axisminmax = [ 0.50, 1.0, -0.05, 1.05 ]  # for consistency

for cmplx_i in range(len(complex_names)):
  plt.plot(pv,papproximation[25,cmplx_i],colors[cmplx_i]+'o-')
  
plt.title('Fraction of p-approximating conformations', weight='bold')  
plt.axis(axisminmax)
plt.xlabel('value of p', weight='bold')
plt.ylabel('fraction well-formed', weight='bold')
plt.legend(short_complex_names,loc='lower left')
plt.savefig("figT25_{}_p_approx.pdf".format(MODE), bbox_inches='tight')
plt.close()

### look at the histograms for completion times for each reaction

filenames=''
i=1
for T in [15, 20, 25, 30, 35, 40, 45, 50, 55, 60]:

  for rxn_i in range(3):  # just the three desired reactions
    plt.figure()
    # histbins[T,rxn_i] = 1000 * histbins[T,rxn_i]   # if ms x scale
    bincenters = [0.5*histbins[T,rxn_i,b]+0.5*histbins[T,rxn_i,b+1] for b in range(30)]
    norm = histdata[T,rxn_i].sum()
    plt.hist(bincenters,weights=histdata[T,rxn_i]/norm,bins=histbins[T,rxn_i])
    axes = plt.gca()
    axes.set_ylim([0,0.20])  # might cut off some
    plt.title('T = {} : {}'.format(T,short_reaction_names[rxn_i]), weight='bold')
    plt.xlabel('reaction time, s', weight='bold')
    #  plt.xlabel('reaction time, ms')
    plt.ylabel('frequency', weight='bold')
    plt.savefig("plot{}.pdf".format(i), bbox_inches='tight')
    plt.close()
    if rxn_i==0 and T in [25,60] and MODE=='publication':
      os.system('cp plot{}.pdf figT{}_publication_hist.pdf'.format(i,T))
    filenames=filenames+' plot{}.pdf'.format(i)
    i=i+1

cmd="pdfjam"+filenames+" --nup 3x10 --outfile fig7_"+MODE+"_hist.pdf --papersize '{6in,15in}'"
os.system(cmd)
cmd="rm"+filenames
os.system(cmd)


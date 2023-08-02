# Build CRNs for Figure 7 (Entropy-driven catalyst, Zhang et al, Science 2007)
# -- Usage is parallel to fig7_plotter.py, but automatically reads all temperatures, and requires a directory name.

# Usage:
#   python fig7_build_crns.py . demo
#   python fig7_build_crns.py publication_data randomT

# Note: Uses the *.csv files produced by fig7_analyze.py

import sys, os, csv
import numpy as np

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

species=['Catalyst','Substrate','Fuel','Intermediate','Output','Signal','Waste'] 

for temp in range(0,100):
  ANALYSIS_PATH = '{}fig7_T{:d}_{}_analysis.csv'.format(DIR,temp,MODE)
  if os.path.exists(ANALYSIS_PATH):
    print("Loading data for temperature T = {:d} and {} mode.".format(temp,MODE))
    temps.append(temp)

    # produce input file for Stefan Badelt's pilsimulator in Python

    CRN_PATH = 'fig7_T{:d}_{}.crn'.format(temp,MODE)
    out = open(CRN_PATH, 'w')
    out.write('# CRN for Zhang 2007 ({}) at T = {} including all observed reactions from simulations\n\n'.format(MODE,temp))

    i=0
    with open(ANALYSIS_PATH, 'r') as f:
      reader = csv.reader(f)
      for row in reader:
        if len(row)==8 and not row[0].find('->') == -1:
          [reactant_string,product_string]=row[0].split('->')
          reactants = [ [s for s in species if s in c] for c in reactant_string.split('+')]
          reactants = [ r for rs in reactants for r in rs ]
          products  = [ [s for s in species if s in c] for c in product_string.split('+')]
          products  = [ r for rs in products for r in rs ]
          if product_string in [" [Complex(cpx_F:LB:OB)] + [Complex(cpx_SB)]", " [Complex(cpx_F:SB:LB)] + [Complex(cpx_OB)]"]:
            products = ["Leak"]
          if len(products)==0:
            break
          [k1,k1_err,k2,k2_err] = [ float(row[1]), float(row[2]), float(row[3]), float(row[3]) ]
          out.write('reaction [k1 = {} +/- {} /M/s] {} -> I{}\n'.format(k1,k1_err,' + '.join(reactants),i))
          out.write('reaction [k2 = {} +/- {} /s  ] I{} -> {}\n'.format(k2,k2_err,i,' + '.join(products)))
          i=i+1
    out.write('\n')
    out.close()

    # produce input file for David Soloveichik's CRN Simulator in Mathematica

    MMA_PATH = 'fig7_T{:d}_{}_crn.m'.format(temp,MODE)
    out = open(MMA_PATH, 'w')
    out.write('(* CRN for Zhang 2007 ({}) at T = {} including all observed reactions from simulations *)\n\n'.format(MODE,temp))

    out.write('x = 10 * 10^-9\n')
    out.write('KindaCRN = {\n')

    i=0
    with open(ANALYSIS_PATH, 'r') as f:
      reader = csv.reader(f)
      for row in reader:
        if len(row)==8 and not row[0].find('->') == -1:
          [reactant_string,product_string]=row[0].split('->')
          reactants = [ [s for s in species if s in c] for c in reactant_string.split('+')]
          reactants = [ r for rs in reactants for r in rs ]
          products  = [ [s for s in species if s in c] for c in product_string.split('+')]
          products  = [ r for rs in products for r in rs ]
          if product_string in [" [Complex(cpx_F:LB:OB)] + [Complex(cpx_SB)]", " [Complex(cpx_F:SB:LB)] + [Complex(cpx_OB)]"]:
            products = ["Leak"]
          if len(products)==0:
            break
          [k1,k1_err,k2,k2_err] = [ float(row[1]), float(row[2]), float(row[3]), float(row[3]) ]
          out.write('  rxn[{},I{},{}], rxn[I{},{},{}],\n'.format('+'.join(reactants),i,k1,i,'+'.join(products),k2))
          i=i+1

    out.write('  conc[Fuel, 1.3 x], conc[Substrate, 1.0 x]\n')
    out.write('}\n')

    out.write('\n')
    out.close()


print()

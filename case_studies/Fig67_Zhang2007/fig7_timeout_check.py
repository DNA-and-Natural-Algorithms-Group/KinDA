# Looks for and reports which reactions had timeouts.  
# If there were timeouts, then statistics may be influenced, and simulations should be redone with a larger timeout limit if feasible.
# -- Usage is parallel to fig7_analyze.py, but automatically reads all temperatures, and requires a directory name.

# Usage:
#   python fig7_timeout_check.py . demo
#   python fig7_timeout_check.py publication_data publication

# Note: Uses the *.csv files produced by fig7_analyze.py

import sys, os, csv

if not len(sys.argv) == 3: 
  sys.exit('ERROR: Needs 2 arguments: a directory name, and simulation type.')

DIR = sys.argv[1]
MODE = sys.argv[2]

MODELIST = [ 'demo', 'trial', 'publication', 'randomT', 'randomR' ]
if not MODE in MODELIST:
  sys.exit('Unrecognized mode.')

if not DIR[-1] == '/':
  DIR = DIR + '/'

reaction_names=[
  '[Complex(Catalyst)] + [Complex(Substrate)] -> [Complex(Intermediate)] + [Complex(Signal)]',
  '[Complex(Intermediate)] + [Complex(Signal)] -> [Complex(Catalyst)] + [Complex(Substrate)]',
  '[Complex(Fuel)] + [Complex(Intermediate)] -> [Complex(Catalyst)] + [Complex(Output)] + [Complex(Waste)]',
  '[Complex(Fuel)] + [Complex(Substrate)] -> [Complex(cpx_F:LB:OB)] + [Complex(cpx_SB)]',
  '[Complex(Fuel)] + [Complex(Substrate)] -> [Complex(cpx_F:SB:LB)] + [Complex(cpx_OB)]'
  ]
short_reaction_names=[
  'Cat+Sub->Int+Sig',
  'Int+Sig->Cat+Sub',
  'Fuel+Int->Cat+Out+Waste',
  'Fuel+Sub->LeakInt1+Sub',
  'Fuel+Sub->LeakInt2+Out'
  ]

for temp in range(0,100):
  ANALYSIS_PATH = '{}fig7_T{:d}_{}_analysis.csv'.format(DIR,temp,MODE)
  if os.path.exists(ANALYSIS_PATH):
    print "Loading data for temperature T = {:d} and {} mode.".format(temp,MODE)

    with open(ANALYSIS_PATH, 'r') as f:
      reader = csv.reader(f)
      for row in reader:
        if len(row)>0 and not row[0].find('HISTOGRAM') == -1:
          break
        if len(row)>0 and not row[0].find('->') == -1 and row[0] in reaction_names:
          rxn_i = reaction_names.index(row[0])
          successes = int(row[5])
          timeouts  = int(row[6])
          total     = int(row[7])
          if timeouts>0:
            print "T = {}, {}: {} successes, {} timeouts, {} total".format(temp,short_reaction_names[rxn_i],successes,timeouts,total)
print


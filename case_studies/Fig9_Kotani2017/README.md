# KinDA case study: Kotani 2017 - Figure 2

## Configuration note
This case study also illustrates the commandline interface to KinDA via the
script `src/kinda/scripts/KinDA.py`, which is installed together with the KinDA
package. Make sure that the script is accessible to your shell environment by
running:

```sh
KinDA --version
```

## Full system analysis
It is recommended to enumerate separately with `peppercorn` first, to have
explicit control over the enumeration semantics. Note, the output file for the
command below (`kotani2017_F2_enum.pil`) can already be found in this
directory.

```sh
~$ cat kotani2017_F2.pil | peppercorn \
        -c -d --max-complex-size 7 > kotani2017_F2_enum.pil
```

Then, let the commandline executable `KinDA` do the rest (eventually -- this
could take days):

```sh
~$ cat kotani2017_F2_enum.pil | KinDA --verbose \
        --temperature 55 --multistrand-timeout 10 \
        --prob-error-goal 0.1 --rate-error-goal 0.4 \
        --backup kotani2017_F2_ms10_T55.db
```

Take a look at the default parameters (`KinDA --help`). Your initial run will
likely not result in the requested relative error goal for each complex
probability (0.1) and reaction rate (0.4). In subsequent runs, use the following
command to continue your analysis using the previously collected data:

```sh
~$ KinDA -v -T 55 --multistrand-timeout 10 -d kotani2017_F2_ms10_T55.db [options]
```

Note that it is not necessary to specify non-default system parameters of the
imported database (here: `-T 55 --multistrand-timeout 10`), but it is good
practice and will avoid warnings.

## Plot your data
Once your error goals are satisfied, use the output of KinDA for further
analysis:

```sh
~$ KinDA -s 0 -r kotani2017_F2_ms10_T55.db > kotani2017_F2_ms10_T55_kinda.pil
```

For example, you can use the `pilsimulator`, which is an executable provided by
the `peppercornenumerator` library to simulate it: 

```sh
~$ cat kotani2017_F2_ms10_T55_kinda.pil | pilsimulator \
        --atol 1e-13 --rtol 1e-13 --mxstep 100000 \
        --t0 0.01 --t-log 10000 --t8 360000 \
        --p0 S1=10 S2=10 R=20 C1=1 \
        --labels C1 S1 S2 R D I1 e51 P1 P2 P3 \
        --pyplot kotani2017_F2_ms10_T55_kinda.pdf
```


## Tips & Tricks
### One count-by-complex reaction
The previous example uses the default `ordered-complex` macrostate definition,
but that does not make much sense for the reaction {S2 + R -> e51}. To collect
data on a single reaction, remove all other condensed reactions and resting
macrostates from the input file (e.g. as shown in
`kotani2017_F2_rxn_RpS2_to_e51.pil`).

```sh
~$ cat kotani2017_F2_rxn_RpS2_to_e51.pil | KinDA -v \
        -T 55 --multistrand-timeout 10 --macrostate-mode count-by-domain \
        -b kotani2017_F2_rxn_RpS2_to_e51_T55.db \
        -e 0.4 --prob-max-sims 0
```


### A rarely successful reaction 
The reaction {I1 + P1 -> S1 + C1} is rarely successful. To speed up your data
collection, you can run multiple invocations in parallel:

```sh
~$ cat kotani2017_F2_slow_rxn.pil | KinDA -v \
        -T 55 --multistrand-timeout 10 -b kotani2017_F2_ms10_T55_1.db \
        --prob-max-sims 0 --rate-max-sims 1000000 --rate-batch-size ...
~$ cat kotani2017_F2_slow_rxn.pil | KinDA -v \
        -T 55 --multistrand-timeout 10 -b kotani2017_F2_ms10_T55_2.db \
        --prob-max-sims 0 --rate-max-sims 1000000 --rate-batch-size ...
```

then merge the files with the original database:

```sh
~$ KinDA -v -s 0 \
    -T 55 --multistrand-timeout 10 \ # to avoid warnings
    -d kotani2017_F2_ms10_T55.db \
    --merge kotani2017_F2_ms10_T55_1.db kotani2017_F2_ms10_T55_2.db ...
```

Beware that the merge option loads data from KinDA databases into your present
system. It will (for good reason) not complain if you merge KinDA files that
have conflicting system parameters!!! You(!) have to make sure that the merge
operation does not lead to skewed data.

### Simulating a reaction with too few successful runs
When KinDA does not have enough simulation data (specifically, 0 or 1 successful
runs), then it will report an infinite error bar on some reaction rates. In the
case of no successful runs, the KinDA script will output a PIL file that
comments out the insufficiently characterized reaction; however, in the case of
a single successful run, the reaction is reported despite the infinite error bar
for k2. While this is reasonable, the released version of `pilsimulator` will
crash because it cannot parse the error bar value. If you encounter this error,
you can either run more simulations, or edit the PIL file by hand to replace
`inf` by a numerical value for the offending reaction. This bug will be fixed in
future versions.

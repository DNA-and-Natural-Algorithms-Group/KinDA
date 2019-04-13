# KinDA case study: Kotani 2017 - Figure 2

## Configuration note

This case study also illustrates the commandline interface to KinDA, via the script `scripts/KinDA`,
which must be in your path for the following commands to work without modification.  E.g. in `bash`,

```sh
PATH=~/KinDA/scripts:$PATH
```

## Full system analysis

It is recommended to enumerate separately with `peppercorn` first, to have
explicit control over the enumeration semantics. Note, the output file for the
command below (`kotani2017_F2_enum.pil`) can already be found in this
directory.

```sh
~$ cat kotani2017_F2.pil | peppercorn -c -d --max-complex-size 7 > kotani2017_F2_enum.pil
```

Then, let the commandline executable `KinDA` do the rest (eventually -- this could take days):

```sh
~$ cat kotani2017_F2_enum.pil | KinDA --verbose \
                                      --temperature 55 --multistrand-timeout 10 \
                                      --error-goal 0.4 \
                                      --backup kotani2017_F2_ms10_T55.db 
```

Take a look at the default parameters (`KinDA --help`). Your initial run will
likely not result with the requested relative error goal (0.4) for each complex
probability or reaction rate. In subsequent runs, use the following command to
continue your analysis using the previously collected data: 

```sh
~$ KinDA -v -T 55 --multistrand-timeout 10 -d kotani2017_F2_ms10_T55.db  [options]
```

Note it is not necessary to specify non-default system parameters of the imported
database (here: `-T 55 --multistrand-timeout 10`), but it is good practice and will
avoid warnings.

## Plot your data
Once your error goals are satisfied, use the output of KinDA for further analysis:

```sh
~$ KinDA -s 0 -r kotani2017_F2_ms10_T55.db > kotani2017_F2_ms10_T55_kinda.pil
```

For example, you can use the `pilsimulator`, which is an executable provided by
the `peppercornenumerator` library to simulate it: 

```sh
~$ cat kotani2017_F2_ms10_T55_kinda.pil | pilsimulator --atol 1e-13 --rtol 1e-13 --mxstep 100000 \
                                                       --t0 0.01 --t-log 10000 --t8 360000 \
                                                       --pyplot-labels C1 S1 S2 R D I1 e51 P1 P2 P3 \
                                                       --p0 S1=10e-9 S2=10e-9 R=20e-9 C1=1e-9 \
                                                       --pyplot kotani2017_F2_ms10_T55_kinda.pdf
```


## Tips & Tricks
### One count-by-complex reaction
The previous example uses the default `ordered-complex` macrostate definition,
but that does not make much sense for the reaction {R + S2 -> e51}. To collect
data on a single reaction, remove all other condensed rections and resting
macrostates from the input file (e.g. as shown in
`kotani2017_F2_rxn_RpS2_to_e51.pil`). 

```sh
~$ cat kotani2017_F2_rxn_RpS2_to_e51.pil | KinDA -v -T 55 --multistrand-timeout 10 \
                                              --macrostate-mode count-by-domain \
                                              -b kotani2017_F2_rxn_RpS2_to_e51_T55.db \
                                              -e 0.4 \
                                              --prob-max-sims 0 
```


### A rarely successful reaction 
The reaction {I1 + P1 -> S1 + C1} is rarely successful. To speed up your data
collection, you can run multiple invocations in parallel:

```sh
~$ cat kotani2017_F2_slow_rxn.pil | KinDA -v -T 55 --multistrand-timeout 10 \
                                       -b kotani2017_F2_ms10_T55_1.db \
                                       --prob-max-sims 0 \
                                       --rate-max-sims 1000000 --rate-batch-size ...
~$ cat kotani2017_F2_slow_rxn.pil | KinDA -v -T 55 --multistrand-timeout 10 \
                                       -b kotani2017_F2_ms10_T55_2.db \
                                       --prob-max-sims 0 \
                                       --rate-max-sims 1000000 --rate-batch-size ...
```

then merge the files with the original database:

```sh
~$ KinDA -v -s 0 \
    -T 55 --multistrand-timeout 10 \ # to avoid warnings
    -d kotani2017_F2_ms10_T55.db \
    --merge kotani2017_F2_ms10_T55_1.db kotani2017_F2_ms10_T55_2.db ...
```

Beware that the merge option loads data from kinda databases into your present
system. It will (for good reason) not complain if you merge kinda files that
have conflicting system parameters!!! You(!) have to make sure that the merge
operation does not lead to skewed data.



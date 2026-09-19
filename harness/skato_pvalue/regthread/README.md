# Region-level parallelism benchmark

Answers: for BURDEN, does splitting the group file across PROCESSES beat
threading inside one process, and how close does either get to linear?

## Why this needs an idle machine

The first run of this measurement (2026-09-19) was taken on a 4-core laptop at
load ~3 with an editor and browser running. Measured cores actually obtained:

| arm | wall | user | cores busy |
|---|---|---|---|
| serial | 30.27 s | 25.43 s | **0.84** |
| 4 threads | 19.74 s | 37.43 s | **1.89** |
| 4 processes | 14.98 s | 48.60 s | **3.24** |

No arm got four cores; the serial baseline did not get one. Those wall times
measure how aggressively each arm grabs CONTENDED cpu, not how it scales, and
they named the wrong winner. **User time is the load-insensitive number** --
25.43 / 37.43 / 48.60 s says threads inflate total CPU 1.47x and processes
1.91x, which on an idle 4-core box predicts threads ~2.7x and processes ~2.1x,
i.e. the opposite ranking. `bench.sh` therefore REFUSES to run above 1-min load
0.5 or below 4 physical cores, and reports user time first.

## Fixture (~28 MB)

Build it once from the sim50k fixture, which is NOT in this repo (it lives in
the separate SAIGE-sim tree under a gitignored `harness/work/`), then copy the
28 MB result to the target machine:

```sh
D=<path-to>/SAIGE-sim/harness/work/sim50k
B=harness/work/benchkit && mkdir -p $B
awk '$2=="var"{for(i=3;i<=NF;i++) print $i}' $D/groupFile.txt | sort -u > $B/keep.txt
.pixi/envs/default/bin/plink2 --pfile $D/chr20_sim_pgen --extract $B/keep.txt \
    --make-pgen --out $B/geno
cp $D/step1_y_binary.rda $D/step1_y_binary.varianceRatio.txt $D/groupFile.txt \
   $D/pedigree_truth.sparseGRM.mtx $D/pedigree_truth.sparseGRM.mtx.sampleIDs.txt $B/
```

Subsetting the pgen to the 37,765 markers the group file names drops it from
389 MB to 10 MB and is inert: **verified byte-identical** to the full-fixture
run (3,345 rows, `cmp` clean). The remaining bulk is `step1_y_binary.rda` at
12 MB.

## Run

```sh
pixi install && bash harness/build.sh      # once per machine
bash harness/skato_pvalue/regthread/bench.sh <fixture-dir> [reps]
```

Each rep also checks that the four quarters reconstruct the whole-file result
and says so; if the row counts are zero or unequal it says the timings are not
valid. That guard exists because the first version of this benchmark failed
outright (missing genotype arguments), produced no output, and its merge check
reported "IDENTICAL" while comparing two EMPTY files.

## What to report back

The three user-time figures and the three wall figures, plus the cores-busy
ratio (user/wall) for each arm. If cores-busy for the 4-way arms is well below
4, the box was not idle and the wall numbers should be discarded.

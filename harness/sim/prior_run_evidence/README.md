# Evidence from the prior (stalled) N=10,000 chr20 run

Logs and metadata only -- the 6.5 GB of generated data was reclaimed on
2026-09-01 to make room for the N=50,000 fixture. Everything here is
regenerable from ../01_simulate.py and ../02_build_dataset.py at seed 7.

Kept because these are the record of WHY that run stalled, and the numbers are
load-bearing for sizing the replacement:

- `01_simulate.log`  -- 1,814,359 trees, 2,972,707 sites, 64.44 Mb, 18,000
                        haplotypes, ~24 min. Gives 46.1k sites/Mb raw.
- `02_build.log`     -- 2,267,192 kept variants, 49,318 coding (2.2%) across
                        550 genes = 89.7 coding variants/gene. 35.2k/Mb after
                        the biallelic + polymorphic-in-EUR filter.
- `grm_check_*.txt`  -- THE reason it stalled: at relatednessCutoff 0.05 the
                        genotype-estimated sparse GRM found 1,115,586 pairs
                        against 2,493 true, i.e. 83x too dense, and the cutoff
                        search selected 0.4. Not a code bug -- the Tennessen
                        model's bottleneck-then-growth gives a large sample
                        genuine background relatedness, and ~2,700 effectively
                        independent markers put a 0.05 cutoff only ~2 SD out.
- `sim_meta.json`    -- the realised parameters of that run.

The replacement fixture avoids this by writing the sparse GRM analytically from
the known pedigree, and treating createSparseGRM's recovery of that truth as a
measured control rather than a prerequisite.

# --nThreads marker scan vs the serial goldens: the two heavy workloads.
#
#   P7  multi-chunk, 5858 markers (chunking happens inside each shard)
#   P11 the AF>0.5 path, post-Stage-B
#
# plus the repeat-identity check (P1, P7 at n=4). Split out of
# test-step2-threads.R purely to shorten the parallel critical path -- see the
# byte-exactness rationale in that file's header, which applies unchanged.
# Scaffolding is in helper-threads.R.
#
# Tier: standard (runs the CLI).

for (key in c("P7", "P11")) {
  local({
    k <- key
    test_that(paste0("threaded marker scan reproduces golden ", k), {
      skip_unless_level("standard")
      skip_if_no_goldens()
      golden <- file.path(saige_goldens, "step2", th_workloads[[k]]$golden)
      for (nt in TH_NTHREADS) {
        out <- th_run(k, nt)
        expect_output_matches_golden(golden, out)
      }
    })
  })
}

test_that("repeated threaded runs are mutually identical (P1, P7 at n=4)", {
  skip_unless_level("standard")
  skip_if_no_goldens()
  # Fork scheduling is nondeterministic; if it could reach the output bytes at
  # all, three repeats of the same command would eventually disagree.
  for (k in c("P1", "P7")) {
    outs <- vapply(1:3, function(r) th_run(k, 4, tag = paste0("_rep", r)),
                   character(1))
    expect_true(files_byte_identical(outs[1], outs[2]))
    expect_true(files_byte_identical(outs[1], outs[3]))
  }
})

# Self-test of harness/timeout.sh, the per-arm watchdog harness/tsan.sh and
# harness/guardmalloc.sh run every arm under. It is here rather than described
# because its whole job is a failure mode that is invisible when it goes wrong:
# a watchdog that does not actually kill leaves the harness hanging exactly as
# it hung before, and reports nothing.
#
# The case that motivated it: on 2026-09-07 a planted TSAN arm segfaulted, R's
# segfault handler printed its banner and then blocked on the interactive abort
# prompt, and tsan.sh sat on it for 4h32m. So the assertions are (a) the status
# of a command that exits, cleanly or by a signal, is reported unchanged, (b) a
# command that IGNORES SIGTERM is still killed and reported as 124, (c) so are
# the children it forked -- SkatoParallel's pool outlives a TERM to the parent
# alone -- and (d) stdin is /dev/null, so a prompt gets EOF instead of waiting.
#
# The hang is observed through files the processes touch, not through pgrep: a
# pattern general enough to match `sleep` also matches whatever else the box is
# running, and this suite runs eight processes of its own. Tier: standard.

test_that("run_with_timeout reports statuses and kills what it must", {
  skip_unless_level("standard")
  if (.Platform$OS.type != "unix") skip("bash watchdog")
  helper <- file.path(saige_repo, "harness", "timeout.sh")
  skip_if_not(file.exists(helper), "harness/timeout.sh missing")

  d <- tempfile("tmo_")
  dir.create(d)
  hung <- function(x) file.exists(file.path(d, x))
  script <- tempfile(fileext = ".sh")
  writeLines(c(
    "set -uo pipefail",
    paste0("source ", shQuote(helper)),
    paste0("D=", shQuote(d)),
    "run_with_timeout 10 $D/a bash -c 'exit 0';       echo \"clean=$RWT_STATUS\"",
    "run_with_timeout 10 $D/b bash -c 'exit 3';       echo \"nonzero=$RWT_STATUS\"",
    "run_with_timeout 10 $D/c bash -c 'kill -SEGV $$'; echo \"signal=$RWT_STATUS\"",
    "run_with_timeout 10 $D/d bash -c 'read -r x';    echo \"stdin=$RWT_STATUS\"",
    # the parent ignores SIGTERM; the child would touch child_end 4 s from now
    # `d` is interpolated by paste0, not left as $D: the command below sits
    # inside single quotes, where the shell expands nothing
    paste0("run_with_timeout 2 $D/e bash -c '",
           "bash -c \"touch ", d, "/child_begin; sleep 4; touch ", d,
           "/child_end\" & ",
           "touch ", d, "/parent_begin; trap \"\" TERM; sleep 120'; ",
           "echo \"hang=$RWT_STATUS\""),
    "sleep 5",
    "echo done=yes"
  ), script)
  out <- system2("bash", script, stdout = TRUE, stderr = FALSE)
  got <- strsplit(grep("=", out, value = TRUE), "=", fixed = TRUE)
  got <- setNames(vapply(got, `[`, "", 2L), vapply(got, `[`, "", 1L))

  expect_identical(got[["clean"]], "0")
  expect_identical(got[["nonzero"]], "3")
  expect_identical(got[["signal"]], "139")   # SIGSEGV, reported as 128 + signal
  # `read` returns 1 at end of file; had stdin been a terminal it would have
  # waited and come back 124 instead
  expect_identical(got[["stdin"]], "1")
  expect_identical(got[["hang"]], "124")
  # the controls: both processes really did start, so their absence five
  # seconds later is the watchdog and not a command that never ran
  expect_true(hung("parent_begin"))
  expect_true(hung("child_begin"))
  expect_false(hung("child_end"))
})

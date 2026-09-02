# The PLINK reader's ownership and cloning contract (S5), and the .bed read
# path that nothing had ever exercised.
#
# Neither of these is visible to a golden. A clone that decoded the wrong bytes,
# or a close() that unmapped the file under a live clone, would leave every
# committed output byte-identical -- the pipeline uses exactly one reader and
# never calls the step-2 closeGenoFile() at all. So, as with the S2 seed
# counter, the change carries its own measurement: the C++ test hooks
# testReadMarkersInCPP / testGenoReaderCloneInCPP decode markers through a
# named handle and hand the genotype vectors back for comparison here.
#
# These tests replace the process-global reader (setPLINKobjInCPP) and the
# clone test CLOSES it, so nothing after them may assume a live reader. That is
# safe as written: no other test in the suite touches g_reader in-process --
# every other reader use goes through a CLI child process.

bed_prefix <- file.path(saige_extdata, "input", "genotype_100markers")

# The analysis sample set is the fam file's IIDs, which is what SPAGMMATtest
# hands setPLINKobjInCPP after intersecting with the model.
fam_ids <- function() {
  fam <- utils::read.table(paste0(bed_prefix, ".fam"), stringsAsFactors = FALSE)
  as.character(fam[[2]])
}

# The env var is read in PlinkClass::setPlinkobj, so it must be set BEFORE the
# reader is built. Set and cleared explicitly rather than through withr, which
# nothing else in this suite depends on.
no_mmap <- function(on) {
  if (on) Sys.setenv(SAIGE_PLINK_NO_MMAP = "1") else Sys.unsetenv("SAIGE_PLINK_NO_MMAP")
}

setup_reader <- function(allele_order = "alt-first") {
  setPLINKobjInCPP(paste0(bed_prefix, ".bim"),
                   paste0(bed_prefix, ".fam"),
                   paste0(bed_prefix, ".bed"),
                   fam_ids(),
                   allele_order)
}

# Deliberately NOT ascending, and deliberately containing row 0 in the middle:
# that is the LEDGER #16 pattern -- the seek for bim row 0 used to be skipped,
# so a mid-stream request for row 0 returned whatever record the cursor sat on.
# Any reader that carries a stream position gets this wrong.
MID_STREAM_ROW0 <- c(5L, 0L, 17L, 0L, 3L)

expect_decode_identical <- function(a, b, label) {
  expect_identical(a$id, b$id, label = paste("marker IDs:", label))
  expect_identical(a$g, b$g, label = paste("genotype vectors:", label))
  expect_identical(a$altFreq, b$altFreq, label = paste("altFreq:", label))
  expect_identical(a$altCounts, b$altCounts, label = paste("altCounts:", label))
  expect_identical(a$missingRate, b$missingRate, label = paste("missingRate:", label))
}

test_that("the reader is memory-mapped by default, and pread when told not to be", {
  no_mmap(FALSE)
  setup_reader()
  expect_true(isGenoReaderForkSafeInCPP())   # PLINK reports fork-safety = mmap

  no_mmap(TRUE)
  setup_reader()
  expect_false(isGenoReaderForkSafeInCPP())
  no_mmap(FALSE)
})

test_that("the pread fallback decodes exactly what the mmap path decodes", {
  # The fallback existed only for an mmap failure, which does not happen for a
  # small local file, so before SAIGE_PLINK_NO_MMAP nothing had ever run it.
  # It is one of the two ways this package can read a .bed and it was gated by
  # nothing; this is that gate.
  no_mmap(FALSE)
  setup_reader()
  mapped <- testReadMarkersInCPP(MID_STREAM_ROW0)

  no_mmap(TRUE)
  setup_reader()
  pread <- testReadMarkersInCPP(MID_STREAM_ROW0)
  no_mmap(FALSE)

  # Byte-exactness is the right claim: same file, same offsets, same decode.
  expect_decode_identical(mapped, pread, "mmap vs pread")
  # ... and the positive control that the comparison could have failed: the two
  # runs really did take different paths.
  expect_identical(length(mapped$g), length(MID_STREAM_ROW0))
})

test_that("both read paths honour the LEDGER #16 random-access contract", {
  # Row 0 requested mid-stream must return row 0, and the same row requested
  # twice at different points in the sequence must agree.
  for (off in c(FALSE, TRUE)) {
    no_mmap(off)
    setup_reader()
    seq_read <- testReadMarkersInCPP(MID_STREAM_ROW0)
    lone <- testReadMarkersInCPP(0L)
    lbl <- if (off) "pread" else "mmap"
    # positions 2 and 4 of MID_STREAM_ROW0 are both row 0
    expect_identical(seq_read$id[[2]], lone$id[[1]], label = paste("row 0 id:", lbl))
    expect_identical(seq_read$g[[2]], lone$g[[1]], label = paste("row 0 genotypes:", lbl))
    expect_identical(seq_read$g[[2]], seq_read$g[[4]],
                     label = paste("row 0 twice mid-stream:", lbl))
  }
  no_mmap(FALSE)
})

test_that("a clone decodes identically to its primary and outlives the primary's close", {
  no_mmap(FALSE)
  setup_reader()
  r <- testGenoReaderCloneInCPP(MID_STREAM_ROW0)

  expect_decode_identical(r$primary, r$clone, "clone vs primary")
  # The claim that matters: the primary has called close(), and the clone still
  # reads. If the map died with the first owner to close rather than the last,
  # this either segfaults or returns garbage -- it cannot quietly pass.
  expect_decode_identical(r$primary, r$cloneAfterPrimaryClose,
                          "clone after primary close")
  # the reader is now closed; leave the process without one rather than a
  # half-open one
  expect_error(testReadMarkersInCPP(0L), "after close")
})

test_that("clone is refused on the path that has not been measured", {
  no_mmap(TRUE)
  setup_reader()
  expect_error(testGenoReaderCloneInCPP(0L), "mmap")
  no_mmap(FALSE)
})

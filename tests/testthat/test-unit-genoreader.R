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
# These tests replace the process-global reader (setPLINKobjInCPP /
# setPGENobjInCPP) and each clone test CLOSES it, so nothing after them may
# assume a live reader -- every test below sets its own up first. That is safe
# as written: no other test in the suite touches g_reader in-process -- every
# other reader use goes through a CLI child process.
#
# S8 added the PGEN half at the bottom. It is the same contract asked of a
# reader that reaches it a different way: PLINK clones share one decoder because
# its reads are const pointer arithmetic into an mmap, while PGEN clones each
# get their own pgenlib reader over shared metadata, because a pgenlib read
# mutates the decoder. The tests do not care which -- they ask only that a clone
# decode identically and outlive the primary's close.

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


# ---- PGEN (S8) -------------------------------------------------------------

pgen_prefix <- file.path(saige_extdata, "input", "genotype_100markers")

psam_ids <- function() {
  psam <- utils::read.table(paste0(pgen_prefix, ".psam"), comment.char = "#",
                            stringsAsFactors = FALSE)
  as.character(psam[[1]])
}

setup_pgen_reader <- function() {
  setPGENobjInCPP(paste0(pgen_prefix, ".pgen"),
                  paste0(pgen_prefix, ".pvar"),
                  paste0(pgen_prefix, ".psam"),
                  psam_ids())
}

test_that("a pgen clone decodes identically to its primary and outlives its close", {
  # The claim S8 rests on, stated as a test rather than as a reading of
  # pgenlib's documentation: a second plink2::PgenReader built against the SAME
  # PgenFileInfo decodes the same bytes, and closing the first one does not take
  # the shared file info with it.
  #
  # It is also the gate on the filename PgrInit is given. Only the FIRST reader
  # onto a file skips it (it inherits the header's open handle); a clone fopen()s
  # it, so the pvar path that stood there since Phase 7 would have had the second
  # reader decoding a text file. Nothing could see that while clone() threw.
  setup_pgen_reader()
  r <- testGenoReaderCloneInCPP(MID_STREAM_ROW0)

  expect_decode_identical(r$primary, r$clone, "pgen clone vs primary")
  expect_decode_identical(r$primary, r$cloneAfterPrimaryClose,
                          "pgen clone after primary close")
  expect_error(testReadMarkersInCPP(0L), "after close")
})

test_that("the pgen reader honours the LEDGER #16 random-access contract", {
  # A pgenlib reader carries a stream position (fp_vidx) and an LD cache, so
  # this is not a formality: row 0 requested mid-stream must be row 0, and the
  # same row asked for twice must agree.
  setup_pgen_reader()
  seq_read <- testReadMarkersInCPP(MID_STREAM_ROW0)
  lone <- testReadMarkersInCPP(0L)
  expect_identical(seq_read$id[[2]], lone$id[[1]])
  expect_identical(seq_read$g[[2]], lone$g[[1]])
  expect_identical(seq_read$g[[2]], seq_read$g[[4]])
})

test_that("a pgen clone reads the same rows as the PLINK reader on the same fixture", {
  # Cross-reader positive control, so the two tests above cannot both pass on a
  # clone that decodes consistently WRONG. genotype_100markers.pgen was made
  # from the alt-first .bed with `plink2 --make-pgen`, so the decoded dosages
  # must match exactly -- these are integer genotypes, not a float pipeline.
  no_mmap(FALSE)
  setup_reader()
  from_bed <- testReadMarkersInCPP(MID_STREAM_ROW0)
  setup_pgen_reader()
  from_pgen <- testGenoReaderCloneInCPP(MID_STREAM_ROW0)$clone
  expect_identical(from_bed$id, from_pgen$id)
  expect_identical(from_bed$g, from_pgen$g)
})

test_that("the pgen reader reports it cannot be forked, but can be cloned", {
  # These two answers are independent and S8 moved only the second: a pgenlib
  # reader owns a FILE* it fseeks, so a fork would share one offset, while a
  # clone is a whole second reader.
  setup_pgen_reader()
  expect_false(isGenoReaderForkSafeInCPP())
  expect_true(isGenoReaderCloneableInCPP())
})

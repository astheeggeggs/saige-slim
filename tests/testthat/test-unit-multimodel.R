# Unit tests for the multi-model plumbing in R/MultiModel.R: manifest parsing,
# the suffix rule, sample-set grouping and the restart agreement check. Fast
# tier -- none of these runs the CLI.

write_manifest <- function(lines) {
  path <- tempfile(fileext = ".txt")
  writeLines(lines, path)
  path
}

# Two real files to point at, so "the file exists" is not the thing under test.
mm_rda <- file.path(saige_goldens, "step1", "binary_sparseGRM_a.rda")
mm_vr <- file.path(saige_goldens, "step1", "binary_sparseGRM_a.varianceRatio.txt")

test_that("a well-formed manifest parses into names and paths", {
  skip_if_no_goldens()
  path <- write_manifest(c(
    paste("binary", mm_rda, mm_vr),
    paste("quant", mm_rda, mm_vr)
  ))
  m <- SAIGE:::readModelManifest(path)
  expect_identical(nrow(m), 2L)
  expect_identical(m$name, c("binary", "quant"))
  expect_identical(m$modelFile[1], mm_rda)
  expect_identical(m$varianceRatioFile[2], mm_vr)
  # colClasses="character": a phenotype literally called NA or T must survive
  # as a name and not become a logical.
  path2 <- write_manifest(c(paste("NA", mm_rda, mm_vr), paste("T", mm_rda, mm_vr)))
  expect_identical(SAIGE:::readModelManifest(path2)$name, c("NA", "T"))
})

test_that("manifest errors name what is wrong", {
  skip_if_no_goldens()
  expect_error(SAIGE:::readModelManifest(tempfile()), "does not exist")
  expect_error(
    SAIGE:::readModelManifest(write_manifest(paste("a", mm_rda))),
    "must have exactly 3"
  )
  expect_error(
    SAIGE:::readModelManifest(write_manifest(character(0))),
    "has no rows|no lines available"
  )
  expect_error(
    SAIGE:::readModelManifest(write_manifest(c(
      paste("a", mm_rda, mm_vr), paste("a", mm_rda, mm_vr)
    ))),
    "must be unique"
  )
  expect_error(
    SAIGE:::readModelManifest(write_manifest(paste("a/b", mm_rda, mm_vr))),
    "must not contain"
  )
  expect_error(
    SAIGE:::readModelManifest(write_manifest(paste("a", "/nope/x.rda", mm_vr))),
    "model file .* does not exist"
  )
  # The literal NA in column 3 is what step 1 writes when it skipped variance
  # ratio estimation; single-variant tests cannot run without one, and the
  # message has to say so rather than complain about a missing file.
  expect_error(
    SAIGE:::readModelManifest(write_manifest(paste("a", mm_rda, "NA"))),
    "skipped variance-ratio estimation"
  )
})

test_that("the manifest and --GMMATmodelFile are mutually exclusive", {
  skip_if_no_goldens()
  path <- write_manifest(paste("a", mm_rda, mm_vr))
  expect_error(
    SAIGE:::resolveModelSpecs(mm_rda, mm_vr, path),
    "mutually exclusive"
  )
})

test_that("the suffix rule: only a MULTI-row manifest renames the output", {
  skip_if_no_goldens()
  legacy <- SAIGE:::resolveModelSpecs(mm_rda, mm_vr, "")
  expect_identical(length(legacy), 1L)
  expect_identical(legacy[[1]]$suffix, "")

  one <- SAIGE:::resolveModelSpecs("", "", write_manifest(paste("a", mm_rda, mm_vr)))
  expect_identical(length(one), 1L)
  # A one-row manifest is byte- AND name-identical to the legacy run. This is
  # upstream's rule and it is what makes the P1 golden usable as a negative
  # control for the whole feature.
  expect_identical(one[[1]]$suffix, "")

  two <- SAIGE:::resolveModelSpecs("", "", write_manifest(c(
    paste("a", mm_rda, mm_vr), paste("b", mm_rda, mm_vr)
  )))
  expect_identical(vapply(two, function(x) x$suffix, character(1)), c("_a", "_b"))
})

test_that("models are grouped by identical() sample vectors, in manifest order", {
  a <- c("s1", "s2", "s3")
  b <- c("s1", "s2")
  g <- SAIGE:::groupModelsBySampleSet(list(a, b, a, b, a))
  expect_identical(g$groupOf, c(1L, 2L, 1L, 2L, 1L))
  expect_identical(length(g$representatives), 2L)
  expect_identical(g$representatives[[1]], a)

  # identical(), not setequal(): the ORDER of sampleID is the reader's output
  # order, so a permutation is a different group and gets its own reader.
  g2 <- SAIGE:::groupModelsBySampleSet(list(a, rev(a)))
  expect_identical(g2$groupOf, c(1L, 2L))

  # One model is one group.
  expect_identical(SAIGE:::groupModelsBySampleSet(list(a))$groupOf, 1L)
})

test_that("a restart with disagreeing index files stops and names them", {
  ok <- list(list(indexChunk = 4), list(indexChunk = 4))
  expect_identical(SAIGE:::restartChunkForModels(ok, c("a.index", "b.index")), 4L)

  bad <- list(list(indexChunk = 4), list(indexChunk = 2))
  expect_error(
    SAIGE:::restartChunkForModels(bad, c("a.index", "b.index")),
    "a.index -> chunk 4.*b.index -> chunk 2"
  )
})

test_that("P > 1 is refused for the two paths that keep per-run state (regions no longer)", {
  skip_if_no_goldens()
  path <- write_manifest(c(paste("a", mm_rda, mm_vr), paste("b", mm_rda, mm_vr)))
  # Each stop must fire BEFORE anything is pushed into C++, which is why they
  # are testable in-process at all: no reader, no model, no output file has been
  # touched when they raise.
  base <- list(
    bedFile = file.path(saige_extdata, "input", "genotype_100markers.bed"),
    SAIGEOutputFile = tempfile(),
    GMMATmodel_varianceRatio_multiTraits_File = path
  )
  # Regions took a manifest from step 4 on (2026-09-06): the run proceeds to
  # the group-file check, which is what fires here on a file that does not
  # exist. The refusal that remains -- models on DIFFERENT sample sets -- needs
  # the models read and is gated end to end in test-step2-region-multimodel.R.
  expect_error(
    do.call(SAIGE::SPAGMMATtest, c(base, list(groupFile = "nope.txt"))),
    "groupFile nope.txt does not"
  )
  expect_error(
    do.call(SAIGE::SPAGMMATtest, c(base, list(is_rewrite_XnonPAR_forMales = TRUE))),
    "male sample indices"
  )
  expect_error(
    do.call(SAIGE::SPAGMMATtest, c(base, list(subSampleFile = "nope.txt"))),
    "one subsample to one model"
  )
})

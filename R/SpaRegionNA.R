# End-of-run accounting for LEDGER #181 (2026-10-16): the saddlepoint R runs
# itself returning NOTHING.
#
# The region path's burden rescale (get_newPhi_scaleFactor_traitType) and the
# conditioning-marker rescale call the C++ saddlepoint kernel, SPA_pval,
# directly from R, and its result went straight into
# SPA_ER_kernel_related_Phiadj_fast_new's `if (VarQ_2 == 0)`. A NaN from the
# kernel -- reachable since the LEDGER #168 returns hand the kernel's value
# back without a finiteness check -- therefore errored OUTSIDE the per-cell
# kernel catch (LEDGER #57) and killed the chromosome's R process; the wrapper
# delivered the surviving chromosomes and exited 0, with `chrStatus` the only
# record (afr HTN, four jobs, 2026-10-16). This is LEDGER #170's shape again
# and the third time a value outside the catch cost a chromosome.
#
# Now a non-finite burden p is (1) dumped when SAIGE_SPA_DUMP_DIR is set (the
# vector the kernel saw, individual-level: analyse it where the data lives),
# (2) counted and named here, and (3) replaced by the normal-approximation p
# with the convergence flag FALSE, so the cell is tested and the region
# survives. The substitution is DISCLOSED, never silent: one line at the end
# of every region run, and `<SAIGEOutputFile>.spaRegionNA.txt` names each
# cell. SAIGE_SPA_PLANT_NA_R=1 is the test plant (every rescale saddlepoint
# returns NA).
spa_region_na_env <- new.env(parent = emptyenv())
spa_region_na_env$rows <- list()
spa_region_na_env$n_dumped <- 0L

spaRegionNAReset <- function() {
  spa_region_na_env$rows <- list()
  spa_region_na_env$n_dumped <- 0L
  invisible(NULL)
}

#' Write the vector a failed rescale saddlepoint saw, when SAIGE_SPA_DUMP_DIR is set.
#' @return the file written, or "" when dumping is off
#' @noRd
spaRegionNADump <- function(mu, g, q, qinv, pval_noadj, n_dropped) {
  dir <- Sys.getenv("SAIGE_SPA_DUMP_DIR", "")
  if (!nzchar(dir)) {
    return("")
  }
  spa_region_na_env$n_dumped <- spa_region_na_env$n_dumped + 1L
  path <- file.path(dir, sprintf("spaRegionNA_%d_%d.tsv", Sys.getpid(), spa_region_na_env$n_dumped))
  con <- file(path, "w")
  on.exit(close(con))
  fmt <- paste0("# kernel=SPA_pval(region rescale) q=%.17g qinv=%.17g ",
                "pval_noadj_log=%.17g n_dropped=%d n=%d")
  hdr <- sprintf(fmt, q, qinv, pval_noadj, as.integer(n_dropped), length(g))
  cat(hdr, "\n", sep = "", file = con)
  cat("mu\tg\n", file = con)
  utils::write.table(data.frame(mu = mu, g = g), con, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  path
}

spaRegionNARecord <- function(kind, region, group, maxMAF, weight, info) {
  if (is.null(info)) {
    return(invisible(NULL))
  }
  spa_region_na_env$rows[[length(spa_region_na_env$rows) + 1L]] <- data.frame(
    Kind = kind, Region = as.character(region), Group = as.character(group),
    max_MAF = as.character(maxMAF),
    Weight = if (is.null(weight) || length(weight) == 0) NA_character_ else as.character(weight),
    n_carriers = as.integer(info$n_carriers), n_pinned_dropped = as.integer(info$n_dropped),
    q_adj = info$q, log_p_normal_reported = info$pval_noadj,
    kernel_value = as.character(info$value), dump = info$dump,
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

spaRegionNATable <- function() {
  if (length(spa_region_na_env$rows) == 0) {
    return(data.frame(Kind = character(), Region = character(), Group = character(),
                      max_MAF = character(), Weight = character(), n_carriers = integer(),
                      n_pinned_dropped = integer(), q_adj = numeric(),
                      log_p_normal_reported = numeric(), kernel_value = character(),
                      dump = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, spa_region_na_env$rows)
}

#' Print the end-of-run line and, if any rescale saddlepoint returned nothing, the sidecar.
#' @noRd
spaRegionNASummary <- function(OutputFile) {
  tab <- spaRegionNATable()
  n <- nrow(tab)
  if (n == 0) {
    cat("Region saddlepoint NA: 0 burden rescales had a saddlepoint that returned nothing.\n")
    return(invisible(tab))
  }
  path <- paste0(OutputFile, ".spaRegionNA.txt")
  utils::write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
  by_kind <- sort(table(tab$Kind), decreasing = TRUE)
  cat("Region saddlepoint NA: ", n, " burden rescale(s) had a saddlepoint that returned nothing ",
      "(NaN/NA) and reported the NORMAL-approximation p for the rescale instead, with the ",
      "convergence flag FALSE (LEDGER #181); before 2026-10-16 this killed the chromosome. ",
      "Each is named in '", path, "'",
      if (any(nzchar(tab$dump))) " with the vector dumped (SAIGE_SPA_DUMP_DIR)" else "",
      ". By kind:\n", sep = "")
  for (k in seq_along(by_kind)) {
    cat("  ", by_kind[[k]], " x ", names(by_kind)[k], "\n", sep = "")
  }
  invisible(tab)
}

# The SKAT-O region results schema, stated once (LEDGER #116).
#
# region_postprocess (R/SAIGE_SPATest_Region.R) assembles a cell's row column
# by column, so before this file the column set existed only implicitly, in
# the data.frame of the first cell that produced one. Two things need it
# WITHOUT a row to read it from:
#   - the header, written when the run starts so a run in which no cell
#     qualifies leaves a parseable, visibly empty file rather than a 0-byte
#     one at exit 0 (LEDGER #116's symptom, reachable by four routes);
#   - the pooled row of a region with nothing tested, "Region Cauchy NA ...",
#     which the BURDEN path (C++ writer) has always written and the SKAT-O
#     path omitted, so the two paths disagreed about how an empty region is
#     reported.
# The order below is the order region_postprocess builds: the cell columns,
# then the conditional columns, then the count columns. It is asserted against
# the first row actually written (SAIGE.Region), so a drift here cannot ship
# silently: every golden run would stop.

regionResultColnames <- function(traitType, isCondition, isMultiWeight) {
  cols <- c("Region", "Group", "max_MAF")
  if (isMultiWeight) cols <- c(cols, "Weight")
  cols <- c(cols, "Pvalue", "Pvalue_Burden", "Pvalue_SKAT")
  if (isMultiWeight) cols <- c(cols, "Pvalue_ACATV", "Pvalue_ACATO")
  cols <- c(cols, "BETA_Burden", "SE_Burden")
  if (isCondition) {
    cols <- c(cols, "Pvalue_cond", "Pvalue_Burden_cond", "Pvalue_SKAT_cond",
              "BETA_Burden_cond", "SE_Burden_cond")
  }
  cols <- c(cols, "MAC")
  if (traitType == "binary") cols <- c(cols, "MAC_case", "MAC_control")
  c(cols, "Number_rare", "Number_ultra_rare")
}

# The header line fwrite() writes for that schema: tab-separated, one newline.
regionHeaderLine <- function(traitType, isCondition, isMultiWeight) {
  paste0(paste(regionResultColnames(traitType, isCondition, isMultiWeight), collapse = "\t"),
         "\n")
}

# The one row a region with nothing tested contributes: its name, Group
# "Cauchy", every other column NA -- the BURDEN writer's row for the same
# region. The NAs are LOGICAL on purpose: rbind() with a numeric column keeps
# the column numeric and with a character column keeps it character, so a
# chunk that holds this row beside real regions prints those regions' numbers
# exactly as it would without it (a character NA would coerce a numeric
# column and change its printed digits).
regionEmptyPooledRow <- function(regionName, traitType, isCondition, isMultiWeight) {
  cols <- regionResultColnames(traitType, isCondition, isMultiWeight)
  row <- as.list(rep(NA, length(cols)))
  names(row) <- cols
  row$Region <- regionName
  row$Group <- "Cauchy"
  as.data.frame(row, stringsAsFactors = FALSE)
}

# Stop if a data.frame of rows about to be written does not carry the schema
# the run's header promised. NULL (a chunk with no rows) passes.
regionAssertSchema <- function(df, traitType, isCondition, isMultiWeight) {
  if (is.null(df)) return(invisible(TRUE))
  expected <- regionResultColnames(traitType, isCondition, isMultiWeight)
  if (!identical(colnames(df), expected)) {
    stop("region results schema drift (LEDGER #116): the rows carry columns [",
         paste(colnames(df), collapse = ", "), "] but the header written at start was [",
         paste(expected, collapse = ", "), "]\n")
  }
  invisible(TRUE)
}

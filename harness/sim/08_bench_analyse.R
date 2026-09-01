# Decompose each step-2 run into startup + per-marker work, and check the two
# readers agree on the association results.
#
# SAIGE prints a cumulative proc.time() after each chunk. Every chunk READS the
# same number of markers (markers_per_chunk), while the number TESTED varies with
# the MAC filter -- so a reader comparison must regress on markers read. The
# slope is seconds per marker read; the intercept at zero markers is startup,
# which is what has to be removed for this to be a fair fight. R^2 is reported:
# if per-chunk cost is not constant the decomposition is not trustworthy and
# should not be quoted.
suppressWarnings(suppressMessages(library(stats)))
dir <- "harness/work/readerbench"
logs <- sort(Sys.glob(file.path(dir, "*_*.log")))
rows <- list()
for (f in logs) {
  b <- sub("\\.log$", "", basename(f)); arm <- sub("_.*", "", b); rep <- sub(".*_", "", b)
  L <- readLines(f, warn = FALSE)
  cum <- as.numeric(sub("^\\s*[0-9.]+\\s+[0-9.]+\\s+([0-9.]+)\\s*$", "\\1",
                        grep("^\\s*[0-9.]+\\s+[0-9.]+\\s+[0-9.]+\\s*$", L, value = TRUE)))
  # The FINAL chunk is partial (210,515 markers = 21 full chunks of 10,000 plus
  # 515) and prints no "Completed X/Y" line, so timing blocks outnumber chunk
  # lines by one. Regress on the chunk ORDINAL over the FULL chunks only: each
  # reads exactly markers_per_chunk, so the slope is seconds per fixed block of
  # markers READ and the intercept is startup. Including the short final chunk
  # would flatten the slope and inflate the apparent startup.
  full <- grep("Completed ([0-9]+)/\\1 markers in the chunk", L, value = TRUE)
  chunk_size <- if (length(full)) as.numeric(sub("^Completed ([0-9]+)/.*", "\\1", full[1])) else NA
  tested <- sum(as.numeric(sub(" markers were tested.*", "", grep("markers were tested", L, value = TRUE))))
  rl <- grep("^real", L, value = TRUE)
  real <- if (length(rl)) as.numeric(sub("^real\\s+", "", tail(rl, 1))) else NA_real_
  nfull <- length(full)
  if (is.na(real) || is.na(chunk_size) || nfull < 3 || length(cum) < nfull) {
    cat(sprintf("%-5s rep %s: %s (%d timings, %d full chunks) -- SKIPPED\n",
                arm, rep, if (is.na(real)) "unfinished" else "unusable",
                length(cum), nfull)); next
  }
  idx <- seq_len(nfull)
  fit <- lm(cum[idx] ~ idx)
  rows[[length(rows) + 1]] <- data.frame(
    arm = arm, rep = rep, real_s = real, chunks = length(cum),
    markers_tested = tested,
    startup_s = unname(coef(fit)[1]),
    us_per_marker = unname(coef(fit)[2]) / chunk_size * 1e6,
    work_s = tail(cum, 1) - unname(coef(fit)[1]),
    r2 = summary(fit)$r.squared)
}
d <- do.call(rbind, rows)
if (is.null(d)) { cat("no usable runs\n"); quit() }
cat("\n== per run ==\n")
print(format(d, digits = 4), row.names = FALSE)

cat("\n== per arm (mean +/- sd) ==\n")
for (a in unique(d$arm)) {
  s <- d[d$arm == a, ]
  cat(sprintf("%-5s n=%d  real %6.2f +/- %4.2f s | startup %5.2f +/- %4.2f s | work %6.2f +/- %4.2f s | %6.3f us/marker | R2 %.5f\n",
              a, nrow(s), mean(s$real_s), sd(s$real_s), mean(s$startup_s), sd(s$startup_s),
              mean(s$work_s), sd(s$work_s), mean(s$us_per_marker), min(s$r2)))
}
if (length(unique(d$arm)) == 2) {
  b <- d[d$arm == "bed", ]; p <- d[d$arm == "pgen", ]
  cat(sprintf("\nratio bed/pgen   wall %.3fx | startup %.3fx | WORK %.3fx (startup removed)\n",
              mean(b$real_s) / mean(p$real_s), mean(b$startup_s) / mean(p$startup_s),
              mean(b$work_s) / mean(p$work_s)))
  tt <- try(t.test(b$work_s, p$work_s), silent = TRUE)
  if (!inherits(tt, "try-error"))
    cat(sprintf("work_s Welch t-test: diff %.2f s, 95%% CI [%.2f, %.2f], p = %.4g\n",
                diff(rev(tt$estimate)), tt$conf.int[1], tt$conf.int[2], tt$p.value))
}

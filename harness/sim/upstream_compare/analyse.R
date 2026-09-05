#!/usr/bin/env Rscript
# Timing analysis over results.tsv written by run_arm.sh.
#   Rscript analyse.R [results.tsv]
# Per (workload, trait, nthreads, models): mean +- sd of user CPU and wall per build,
# main/strip ratio, and a PAIRED t-test on user_s by rep (arms adjacent within a rep).
# For s2single, startup is separated by regressing cumulative proc.time() on the chunk
# ordinal over FULL chunks (harness/sim/analyse_all.R's estimator); R2 < 0.99 is flagged.
# For step1, the PCG-iteration budget and "fitting the NULL model took" are read from the
# log, both arms, so a speed ratio can be reconciled against the pass count.
args <- commandArgs(TRUE)
f <- if (length(args)) args[1] else "/Users/dpalmer/Repositories/SAIGE/harness/work/upstream_compare/results.tsv"
d <- read.delim(f, stringsAsFactors = FALSE)
d <- d[d$rc == 0, ]
chunkfit <- function(log) {
  L <- readLines(log, warn = FALSE)
  cum <- as.numeric(sub("^\\s*[0-9.]+\\s+[0-9.]+\\s+([0-9.]+)\\s*$", "\\1",
           grep("^\\s*[0-9.]+\\s+[0-9.]+\\s+[0-9.]+\\s*$", L, value = TRUE)))
  full <- grep("Completed ([0-9]+)/\\1 markers in the chunk", L, value = TRUE)
  n <- length(full)
  if (n < 3 || length(cum) < n) return(c(startup = NA, work = NA, r2 = NA))
  cs <- as.numeric(sub("^Completed ([0-9]+)/.*", "\\1", full[1]))
  fit <- lm(cum[seq_len(n)] ~ seq_len(n))
  c(startup = coef(fit)[[1]], work = tail(cum, 1) - coef(fit)[[1]], r2 = summary(fit)$r.squared,
    us_marker = coef(fit)[[2]] / cs * 1e6)
}
step1info <- function(log) {
  L <- readLines(log, warn = FALSE)
  it <- sum(as.numeric(sub(".*getPCG1ofSigmaAndVector\\D*([0-9]+)\\s*$", "\\1",
            grep("iter from getPCG1ofSigmaAndVector", L, value = TRUE))), na.rm = TRUE)
  i <- grep("fitting the NULL model took", L)
  fit_user <- if (length(i)) { v <- suppressWarnings(as.numeric(strsplit(trimws(L[i[1] + 2]), "\\s+")[[1]])); v[1] } else NA
  c(pcg_iters = it, fit_user = fit_user)
}
ex <- t(sapply(seq_len(nrow(d)), function(i) {
  x <- c(startup = NA, work = NA, r2 = NA, us_marker = NA, pcg_iters = NA, fit_user = NA)
  if (d$workload[i] == "s2single") { v <- chunkfit(d$log[i]); x[names(v)] <- v }
  if (d$workload[i] == "step1") { v <- step1info(d$log[i]); x[names(v)] <- v }
  x
}))
d <- cbind(d, ex)
d$cell <- paste(d$workload, d$trait, paste0("nt", d$nthreads), d$models)
for (cell in unique(d$cell)) {
  s <- d[d$cell == cell, ]
  cat(sprintf("\n== %s ==\n", cell))
  for (b in c("main", "strip")) {
    x <- s[s$build == b, ]; if (!nrow(x)) next
    cat(sprintf("  %-5s n=%d  user %7.2f +- %5.2f   real %7.2f +- %5.2f   maxRSS %6.0f MB", b, nrow(x),
        mean(x$user_s), if (nrow(x) > 1) sd(x$user_s) else 0, mean(x$real_s), if (nrow(x) > 1) sd(x$real_s) else 0, mean(x$maxRSS_MB)))
    if (x$workload[1] == "s2single") cat(sprintf("   work %6.2f +- %4.2f  startup %4.2f  %5.1f us/marker  minR2 %.4f",
        mean(x$work), if (nrow(x) > 1) sd(x$work) else 0, mean(x$startup), mean(x$us_marker), min(x$r2)))
    if (x$workload[1] == "step1") cat(sprintf("   PCG iters %s  null-fit user %s", paste(unique(x$pcg_iters), collapse = "/"), paste(round(x$fit_user, 1), collapse = "/")))
    cat("\n")
  }
  m <- s[s$build == "main", ]; t_ <- s[s$build == "strip", ]
  reps <- intersect(m$rep, t_$rep)
  if (length(reps)) {
    mu <- m$user_s[match(reps, m$rep)]; tu <- t_$user_s[match(reps, t_$rep)]
    cat(sprintf("  main/strip user CPU = %.3fx  (paired on %d reps", mean(mu) / mean(tu), length(reps)))
    if (length(reps) >= 2) { tt <- t.test(mu, tu, paired = TRUE); cat(sprintf(", diff %.2f s, 95%% CI [%.2f, %.2f], p = %.3g", tt$estimate, tt$conf.int[1], tt$conf.int[2], tt$p.value)) }
    cat(")\n")
    if (s$workload[1] == "s2single") {
      mw <- m$work[match(reps, m$rep)]; tw <- t_$work[match(reps, t_$rep)]
      cat(sprintf("  main/strip work_s   = %.3fx   (startup removed)\n", mean(mw) / mean(tw)))
    }
  }
}
cat("\nlibraries loaded:\n"); print(unique(d[, c("build", "lib", "so_md5")]), row.names = FALSE)

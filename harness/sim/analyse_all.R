# Four arms: {strip, main} x {bed, pgen}. Same fixture, model, GRM; nThreads=1;
# warm cache. Startup separated by regressing each run's cumulative proc.time()
# on the chunk ordinal over FULL chunks only (the last chunk is partial and
# prints no "Completed" line, so including it flattens the slope).
dirs <- c(strip = "/Users/dpalmer/Repositories/SAIGE-sim/harness/work/readerbench",
          main  = "/Users/dpalmer/Repositories/SAIGE-main/bench")
rows <- list()
for (bld in names(dirs)) {
  for (f in sort(Sys.glob(file.path(dirs[[bld]], "*_*.log")))) {
    bn <- sub("\\.log$", "", basename(f))
    arm <- sub("_.*", "", bn); rep <- sub(".*_", "", bn)
    if (!arm %in% c("bed", "pgen")) next
    L <- readLines(f, warn = FALSE)
    cum <- as.numeric(sub("^\\s*[0-9.]+\\s+[0-9.]+\\s+([0-9.]+)\\s*$", "\\1",
              grep("^\\s*[0-9.]+\\s+[0-9.]+\\s+[0-9.]+\\s*$", L, value = TRUE)))
    full <- grep("Completed ([0-9]+)/\\1 markers in the chunk", L, value = TRUE)
    cs <- if (length(full)) as.numeric(sub("^Completed ([0-9]+)/.*", "\\1", full[1])) else NA
    rl <- grep("^real", L, value = TRUE)
    real <- if (length(rl)) as.numeric(sub("^real\\s+", "", tail(rl, 1))) else NA_real_
    ul <- grep("^user", L, value = TRUE)
    user <- if (length(ul)) as.numeric(sub("^user\\s+", "", tail(ul, 1))) else NA_real_
    tested <- sum(as.numeric(sub(" markers were tested.*", "", grep("markers were tested", L, value = TRUE))))
    n <- length(full)
    if (is.na(real) || is.na(cs) || n < 3 || length(cum) < n) next
    fit <- lm(cum[seq_len(n)] ~ seq_len(n))
    rows[[length(rows) + 1]] <- data.frame(build = bld, arm = arm, rep = rep,
      real_s = real, user_s = user, tested = tested,
      startup_s = coef(fit)[[1]], work_s = tail(cum, 1) - coef(fit)[[1]],
      us_marker = coef(fit)[[2]] / cs * 1e6, r2 = summary(fit)$r.squared)
  }
}
d <- do.call(rbind, rows)
cat("\n== per arm ==\n")
agg <- do.call(rbind, lapply(split(d, list(d$build, d$arm), drop = TRUE), function(s)
  data.frame(build = s$build[1], arm = s$arm[1], n = nrow(s), tested = s$tested[1],
             real = mean(s$real_s), real_sd = sd(s$real_s),
             user = mean(s$user_s), work = mean(s$work_s), work_sd = sd(s$work_s),
             startup = mean(s$startup_s), us_marker = mean(s$us_marker), r2 = min(s$r2))))
agg <- agg[order(agg$build, agg$arm), ]
for (i in seq_len(nrow(agg))) with(agg[i, ], cat(sprintf(
  "%-6s %-5s n=%d  tested %d | real %6.2f±%4.2f | user %6.2f | startup %4.2f | WORK %6.2f±%4.2f | %6.1f us/marker | R2 %.4f\n",
  build, arm, n, tested, real, ifelse(is.na(real_sd), 0, real_sd), user, startup,
  work, ifelse(is.na(work_sd), 0, work_sd), us_marker, r2)))
g <- function(b, a) agg$work[agg$build == b & agg$arm == a]
cat("\n== within build: pgen / bed on WORK (startup removed) ==\n")
for (b in c("main", "strip")) if (length(g(b,"bed")) && length(g(b,"pgen")))
  cat(sprintf("  %-6s pgen/bed = %.3f  -> %s is faster by %.1f%%\n", b,
      g(b,"pgen")/g(b,"bed"),
      ifelse(g(b,"pgen") < g(b,"bed"), "pgen", "bed"),
      abs(1 - g(b,"pgen")/g(b,"bed")) * 100))
cat("\n== across builds: how much did each path improve? ==\n")
for (a in c("bed", "pgen")) if (length(g("main",a)) && length(g("strip",a)))
  cat(sprintf("  %-5s main/strip = %.2fx  (strip faster)\n", a, g("main",a)/g("strip",a)))

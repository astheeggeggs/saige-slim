# Reading and replaying a saddlepoint dump (SAIGE_SPA_DUMP_DIR; src/SAIGE_test.cpp
# spaDumpVector). Shared by test-step2-spa-dump.R and test-step2-spa-pinned.R.
# A dump is two header lines of tab-separated key/value pairs and one row per
# sample: gtilde (the vector the kernel saw), mu, y, nz (1 for a carrier), g
# (the raw vector the test was handed) and dropped (LEDGER #174: 1 for a pinned
# sample left out of the CGF; its gtilde is written as the 0 the kernel saw).

spa_dump_kv <- function(line) {
  f <- strsplit(sub("^# ", "", line), "\t", fixed = TRUE)[[1]]
  stats::setNames(f[seq(2, length(f), 2)], f[seq(1, length(f), 2)])
}

read_spa_dump <- function(f) {
  l <- readLines(f, n = 2, warn = FALSE)
  list(h1 = spa_dump_kv(l[1]), h2 = spa_dump_kv(l[2]),
       d = utils::read.table(f, header = TRUE, skip = 2, sep = "\t"))
}

spa_dump_num <- function(kv, k) as.numeric(kv[[k]])

# the replay: the kernel the second line names, on the vector and the flags.
# The fast kernel's exact part is the carriers that were not dropped.
replay_spa_dump <- function(dm) {
  mu <- dm$d$mu
  g <- dm$d$gtilde
  nz <- dm$d$nz == 1
  kept <- nz & dm$d$dropped == 0
  q <- spa_dump_num(dm$h1, "q")
  qinv <- spa_dump_num(dm$h1, "qinv")
  p_na <- spa_dump_num(dm$h2, "pNA")
  tol <- spa_dump_num(dm$h2, "tol")
  quiet <- function(expr) {
    out <- NULL
    utils::capture.output(out <- expr)   # the exported kernels print debugging lines
    out
  }
  if (identical(dm$h2[["kernel"]], "fast")) {
    quiet(SPA_binary_fast(mu, g, q, qinv, p_na, FALSE, g[!nz], g[kept], mu[!nz], mu[kept],
                          spa_dump_num(dm$h2, "NAmu"), spa_dump_num(dm$h2, "NAsigma"), tol))
  } else {
    quiet(SPA_binary(mu, g, q, qinv, p_na, tol, FALSE))
  }
}

# Null calibration for LEDGER #126: checks P1 to P5 of the prediction file under docs/state.
# Run with CALIB_DIR set to the working directory; it reads runs/<trait>_<arm>_g<gate>.txt there.
# Gate 1 against gate 5 per (trait, arm): the cell counts against the prediction
# table (P1), what gate 5 removes (P2), the ECDF of the newly admitted cells by
# MAC stratum with a Clopper-Pearson interval (P3/P4), the positive control (P5).
S <- Sys.getenv("CALIB_DIR")
if (S == "") stop("set CALIB_DIR")
D <- "/Users/dpalmer/Repositories/SAIGE-sim/harness/work/sim50k"
K <- 0.01
num_cols <- c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT", "BETA_Burden", "SE_Burden", "MAC",
              "MAC_case", "MAC_control", "Number_rare", "Number_ultra_rare")

rd <- function(tag) {
  f <- file.path(S, "runs", paste0(tag, ".txt"))
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  x <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "",
                  comment.char = "", colClasses = "character")
  x$line <- readLines(f)[-1]
  for (cc in intersect(num_cols, names(x))) x[[cc]] <- suppressWarnings(as.numeric(x[[cc]]))
  if (!"Weight" %in% names(x)) x$Weight <- "-"
  maf <- sprintf("%g", suppressWarnings(as.numeric(x$max_MAF)))
  x$key <- paste(x$Region, x$Group, maf, x$Weight, sep = "|")
  x
}

rdside <- function(tag) {
  f <- file.path(S, "runs", paste0(tag, ".txt.pooledTests.txt"))
  if (!file.exists(f)) return(NULL)
  read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "",
             comment.char = "")
}

cp <- function(k, n) {
  if (n == 0) c(NA, NA) else binom.test(k, n)$conf.int
}

ecdf_line <- function(p, label, alphas = c(0.05, 0.01), expect = NULL) {
  n <- length(p)
  na <- sum(is.na(p))
  s <- sprintf("%-34s n=%4d NA=%d p==1:%d p==0:%d", label, n, na,
               sum(p >= 1, na.rm = TRUE), sum(p <= 0, na.rm = TRUE))
  for (a in alphas) {
    k <- sum(p <= a, na.rm = TRUE)
    ci <- cp(k, n - na)
    flag <- if (is.na(ci[1])) "" else if (a < ci[1]) " HIGH" else if (a > ci[2]) " low" else ""
    s <- paste0(s, sprintf(" | F(%.2f)=%.3f [%.3f,%.3f]%s", a, k / (n - na), ci[1], ci[2], flag))
  }
  if (!is.null(expect)) s <- paste0(s, sprintf("  exp~%.3f", expect))
  cat(s, "\n")
}

# one cell per (Region, MAC stratum): the largest support, ties by file order
primary <- function(cells) {
  cells$supp <- cells$Number_rare + cells$Number_ultra_rare
  cells <- cells[order(cells$Region, cells$MAC, -cells$supp), ]
  cells[!duplicated(paste(cells$Region, cells$MAC)), ]
}

read_pred <- function(collapsed, raw) {
  f <- if (file.exists(file.path(S, collapsed))) collapsed else raw
  cat("prediction table:", f, "\n")
  x <- read.table(file.path(S, f), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  names(x)[1] <- "gene"
  x
}
pred_gene <- read_pred("predicted_cell_mac_collapsed.tsv", "predicted_cell_mac.tsv")
pred_w <- read_pred("predicted_w25_collapsed.tsv", "predicted_w25.tsv")
truth <- read.table(file.path(D, "truth_genes.tsv"), header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)
causal <- truth$gene[truth$causal == 1]

# P1: every printed MAC is an integer and equals the predicted one
check_p1 <- function(c1, arm) {
  cat(sprintf("P1 non-integer MAC rows: %d\n", sum(c1$MAC != round(c1$MAC), na.rm = TRUE)))
  pred <- if (arm %in% c("wplain", "wgrid", "wgridr")) pred_w else pred_gene
  if (arm == "flex") {
    pk <- pred[grepl("^flex:", pred$group), ]
    pk$key <- paste(pk$gene, "high;medium", sprintf("%g", 0.001),
                    paste0("SIM_", sub("^flex:", "", pk$group)), sep = "|")
  } else {
    pk <- pred[!grepl("^flex:", pred$group), ]
    pk$key <- paste(pk$gene, gsub(":", ";", pk$group), sprintf("%g", pk$maxMAF), "-", sep = "|")
  }
  m <- match(c1$key, pk$key)
  hit <- !is.na(m)
  bad <- hit & c1$MAC != pk$MAC[m]
  cat(sprintf(paste0("P1 cells matched to prediction: %d of %d; MAC mismatches: %d; ",
                     "predicted non-empty cells absent from output: %d\n"),
              sum(hit), nrow(c1), sum(bad, na.rm = TRUE),
              sum(pk$MAC > 0 & !(pk$key %in% c1$key))))
  if (any(bad, na.rm = TRUE)) {
    print(head(cbind(c1[which(bad), c("key", "MAC")], pred = pk$MAC[m[which(bad)]]), 10))
  }
}

# P2: gate 5 is gate 1 minus the MAC 1..4 rows; shared rows byte-identical;
# the sidecar's refusal counts reconcile with the distinct sub-5 cells
check_p2 <- function(c1, c5, g1, g5, tag5) {
  only5 <- setdiff(c5$key, c1$key)
  only1 <- c1[!(c1$key %in% c5$key), ]
  shared <- intersect(c1$key, c5$key)
  diffshared <- sum(c1$line[match(shared, c1$key)] != c5$line[match(shared, c5$key)])
  cat(sprintf(paste0("P2 rows only at g5: %d | rows only at g1: %d (MAC range %s..%s; ",
                     "all in 1..4: %s) | g1 rows with MAC<5 that survive at g5: %d | ",
                     "shared cell rows differing: %d\n"),
              length(only5), nrow(only1), min(only1$MAC), max(only1$MAC),
              all(only1$MAC >= 1 & only1$MAC <= 4),
              sum(c1$MAC < 5 & c1$key %in% c5$key), diffshared))
  lost <- table(factor(only1$MAC, levels = 1:4))
  cat("P2 rows removed at g5 by MAC stratum:", paste(names(lost), lost, sep = ":"), "\n")
  side <- rdside(tag5)
  if (!is.null(side)) {
    # the sidecar counts EVALUATED cells refused; a cell carried forward from an
    # identical (same markers) cell of the previous max_MAF is neither evaluated
    # nor refused, so count distinct (Region, Group, Weight, p, MAC) among the rows
    dkey <- paste(only1$Region, only1$Group, only1$Weight, only1$Pvalue_Burden, only1$MAC)
    dist1 <- only1[!duplicated(dkey), ]
    ref <- tapply(dist1$MAC, dist1$Region, length)
    sr <- setNames(side$Cells_refused_minGroupMAC, side$Region)
    regs <- union(names(ref), names(sr)[sr > 0])
    a <- ifelse(is.na(ref[regs]), 0, ref[regs])
    b <- ifelse(is.na(sr[regs]), 0, sr[regs])
    rules <- table(side$Rule)
    cat(sprintf(paste0("P2 sidecar: %d regions listed, refusal counts disagreeing with the ",
                       "g1 sub-5 rows: %d; rules: %s\n"),
                nrow(side), sum(a != b), paste(names(rules), rules, sep = "=", collapse = ", ")))
    if (any(a != b)) {
      bb <- regs[a != b]
      cat("   disagreeing:", paste(bb, "g1-distinct-sub5=", a[a != b], "sidecar=", b[a != b]), "\n")
    }
  } else {
    cat("P2 sidecar: none written at g5\n")
  }
  ca1 <- g1[g1$Group == "Cauchy", ]
  ca5 <- g5[g5$Group == "Cauchy", ]
  chg <- ca1$Region[ca1$line != ca5$line[match(ca1$Region, ca5$Region)]]
  lostreg <- unique(only1$Region)
  cat(sprintf(paste0("P2 Cauchy rows changed: %d; of which in a region that lost a cell: %d; ",
                     "regions that lost a cell: %d\n"),
              length(chg), sum(chg %in% lostreg), length(lostreg)))
}

# P3/P4: the ECDF of the newly admitted cells, by stratum
check_null <- function(c1, g1, g5, trait) {
  pcol <- "Pvalue_Burden"
  sub <- c1[c1$MAC >= 1 & c1$MAC <= 4, ]
  pr <- primary(sub)
  cat(sprintf("-- ECDF, newly admitted cells (MAC 1..4), %s, PRIMARY one per (region, stratum):\n",
              pcol))
  for (mm in 1:4) {
    expect <- if (trait == "y_null_bin") 1 - (1 - K)^mm else 0.05
    ecdf_line(pr[[pcol]][pr$MAC == mm], sprintf("MAC=%d", mm), expect = expect)
  }
  ecdf_line(pr[[pcol]], "pooled 1..4")
  if (trait == "y_null_bin") {
    for (mm in 1:4) {
      q <- pr[[pcol]][pr$MAC == mm]
      cat(sprintf("   MAC=%d p quantiles 1,5,10,25,50%%: %s ; share p>0.9: %.2f\n", mm,
                  paste(signif(quantile(q, c(0.01, 0.05, 0.1, 0.25, 0.5)), 3), collapse = " "),
                  mean(q > 0.9)))
    }
  }
  cat("-- descriptive, ALL sub-5 rows:\n")
  for (mm in 1:4) ecdf_line(sub[[pcol]][sub$MAC == mm], sprintf("MAC=%d (all rows)", mm))
  if ("Pvalue" %in% names(c1)) {
    cat("-- SKAT-O path, the cell's Pvalue (SKAT-O) on the same primary cells:\n")
    for (mm in 1:4) ecdf_line(pr$Pvalue[pr$MAC == mm], sprintf("MAC=%d SKAT-O p", mm))
  }
  cat("-- reference: cells at MAC 5..20 (admitted at both gates),",
      "primary one per (region, stratum):\n")
  ref <- primary(c1[c1$MAC >= 5 & c1$MAC <= 20, ])
  ecdf_line(ref[[pcol]], "MAC 5..20")
  cat("-- pooled Cauchy row per region, g1 vs g5 (descriptive; atom at 1 expected):\n")
  ca1 <- g1[g1$Group == "Cauchy", ]
  ca5 <- g5[g5$Group == "Cauchy", ]
  pc <- if ("Pvalue" %in% names(ca1)) "Pvalue" else "Pvalue_Burden"
  ecdf_line(ca1[[pc]], "Cauchy g1")
  ecdf_line(ca5[[pc]], "Cauchy g5")
}

# P5: the positive control on the signal trait
check_p5 <- function(c1, c5, g1, g5) {
  ca1 <- g1[g1$Group == "Cauchy", ]
  ca5 <- g5[g5$Group == "Cauchy", ]
  cm <- match(ca1$Region, ca5$Region)
  pc <- if ("Pvalue" %in% names(ca1)) "Pvalue" else "Pvalue_Burden"
  thr <- 0.05 / 210
  d1 <- ca1$Region[!is.na(ca1[[pc]]) & ca1[[pc]] <= thr]
  d5 <- ca5$Region[!is.na(ca5[[pc]]) & ca5[[pc]] <= thr]
  cat(sprintf(paste0("P5 detections at %.2e: g5 %d causal + %d non-causal; ",
                     "g1 %d causal + %d non-causal; lost at g1: %s; gained at g1: %s\n"),
              thr, sum(d5 %in% causal), sum(!(d5 %in% causal)),
              sum(d1 %in% causal), sum(!(d1 %in% causal)),
              paste(setdiff(d5, d1), collapse = ","), paste(setdiff(d1, d5), collapse = ",")))
  r <- ca1[[pc]] / ca5[[pc]][cm]
  ncell1 <- table(c1$Region)[ca1$Region]
  ncell5 <- table(c5$Region)[ca1$Region]
  ncell5[is.na(ncell5)] <- 0
  bound <- ncell1 / pmax(ncell5, 1)
  # the (n+k)/n bound holds for an added cell with p <= 0.5; a cell at p > 0.5
  # pulls the Cauchy statistic down harder, by ~1% on this fixture
  cat(sprintf(paste0("P5 pooled-p ratio g1/g5: median %.3f, range %.3f..%.3f; ",
                     "regions where ratio exceeds (n1/n5): %d of %d that gained cells\n"),
              median(r, na.rm = TRUE), min(r, na.rm = TRUE), max(r, na.rm = TRUE),
              sum(r > bound + 1e-9 & ncell1 > ncell5, na.rm = TRUE), sum(ncell1 > ncell5)))
  sub <- c1[c1$MAC >= 1 & c1$MAC <= 4, ]
  cat(sprintf("P5 sub-5 cells on y_binary: %d rows, min p %.2e, rows with p<1e-4: %d\n",
              nrow(sub), min(sub$Pvalue_Burden, na.rm = TRUE),
              sum(sub$Pvalue_Burden < 1e-4, na.rm = TRUE)))
  ps <- primary(sub)
  for (mm in 1:4) ecdf_line(ps$Pvalue_Burden[ps$MAC == mm], sprintf("MAC=%d (signal trait)", mm))
}

arms <- list(plain = c("y_null", "y_null_bin", "y_binary"),
             flex = c("y_null", "y_null_bin", "y_binary"),
             grid = c("y_null", "y_null_bin", "y_binary"),
             wplain = c("y_null", "y_null_bin"), wgrid = c("y_null", "y_null_bin"),
             wgridr = c("y_null", "y_null_bin"))
for (arm in names(arms)) {
  for (trait in arms[[arm]]) {
    tag5 <- sprintf("%s_%s_g5", trait, arm)
    g1 <- rd(sprintf("%s_%s_g1", trait, arm))
    g5 <- rd(tag5)
    cat("\n==================", trait, arm, "==================\n")
    if (is.null(g1)) {
      cat("MISSING g1 RUN\n")
      next
    }
    nog5 <- is.null(g5)
    if (nog5) {
      cat("NO g5 RUN (the SKAT-O path dies at gate 5 on an all-refused region); P2 skipped,",
          "g5 := g1\n")
      g5 <- g1
    }
    c1 <- g1[g1$Group != "Cauchy", ]
    c5 <- g5[g5$Group != "Cauchy", ]
    cat(sprintf("rows: g1 %d cells + %d Cauchy; g5 %d cells + %d Cauchy\n", nrow(c1),
                sum(g1$Group == "Cauchy"), nrow(c5), sum(g5$Group == "Cauchy")))
    check_p1(c1, arm)
    check_p2(c1, c5, g1, g5, tag5)
    if (trait != "y_binary") check_null(c1, g1, g5, trait) else check_p5(c1, c5, g1, g5)
  }
}

a <- commandArgs(trailingOnly=TRUE); dyn.load(file.path(a[2],"exact.so"))
v <- .Call("compare_exhaustive", as.integer(a[1]))
cat(sprintf("pairs compared: %.0f\nbit-differences: %.0f\n", v[1], v[2]))
if (v[2] > 0) cat(sprintf("first at n=%.0f k=%.0f: Rf_lchoose=%.17g saige=%.17g\n", v[3],v[4],v[5],v[6]))

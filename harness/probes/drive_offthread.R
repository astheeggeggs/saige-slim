a <- commandArgs(trailingOnly=TRUE); dyn.load(file.path(a[2],"exact.so"))
v <- .Call("offthread_check", as.integer(a[1]))
cat(sprintf("worker computed %.0f values, bit-differences vs main thread: %.0f\n", v[1], v[2]))

a <- commandArgs(trailingOnly = TRUE)
dyn.load(file.path(a[3], "loc.so"))
v <- .Call(a[1], as.logical(a[2]))
cat(sprintf("RESULT %s main=%.6f worker=%.6f identical=%s\n",
            a[1], v[1], v[2], identical(v[1], v[2])))

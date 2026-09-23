# SAIGE (strip branch) -- a minimal runtime image.
#
# Two stages. `build` installs the full pixi environment, compiles pgenlib
# and the package. `runtime` starts from a plain base and receives ONLY the
# runtime pixi environment (R and the packages SAIGE imports; no compilers,
# no devtools/testthat, no plink2), the installed package, and the three CLI
# scripts on PATH. No example data, no goldens, no tests: mount your own.
#
#   docker build -f docker/Dockerfile -t saige:X.Y.Z .
#   docker run --rm -v /data:/data -w /data saige:X.Y.Z \
#       step2_SPAtests.R --bedFile=... --SAIGEOutputFile=/data/out.txt
#
# The release workflow (.github/workflows/release.yml) builds this on a
# vX.Y.Z tag and checks the tag against DESCRIPTION's Version.

ARG PIXI_VERSION=0.51.0

# ---------------------------------------------------------------- build
FROM ubuntu:24.04 AS build
ARG PIXI_VERSION
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8 PIXI_NO_PATH_UPDATE=1
RUN apt-get update -qq && apt-get install -y --no-install-recommends curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION=v${PIXI_VERSION} sh \
    && mv /root/.pixi/bin/pixi /usr/local/bin/pixi
WORKDIR /app
# manifest first, so the environment layer is cached across source changes
COPY pixi.toml pixi.lock* ./
RUN pixi install -e default && pixi install -e runtime && rm -rf /root/.cache
COPY DESCRIPTION NAMESPACE ./
COPY src ./src
COPY R ./R
COPY harness/build.sh ./harness/build.sh
COPY extdata/createSparseGRM.R extdata/step1_fitNULLGLMM.R extdata/step2_SPAtests.R extdata/extractNglmm.R ./extdata/
# pgenlib (plink-ng) + the package, installed into the RUNTIME environment's
# R library so the runtime stage can take that tree whole.
RUN SAIGE_WARNINGS_RATCHET=0 SAIGE_R_LIBRARY="$(pwd)/.pixi/envs/runtime/lib/R/library" bash harness/build.sh
# Prune the RUNTIME tree HERE, before it is copied: a layer never shrinks,
# so an rm in the runtime stage would leave the full 1.36 GB underneath.
# conda-forge's r-base on Linux depends on the GCC toolchain (so that
# install.packages can compile), which the image does not do. Drop the
# compilers, sanitizer runtimes, static archives and docs; the runtime
# libraries R links (libgfortran, libstdc++, libgomp) stay. The --help calls
# in the runtime stage and the release workflow's step-2 run are the check.
# The second rm is ~37 MB more, all of it inert at runtime: R's message
# catalogues (LC_ALL=C.UTF-8 here, so English either way), pixi's conda-meta
# bookkeeping (no pixi binary in the runtime stage), and the Rcpp/RcppEigen/
# RcppArmadillo headers, which exist for COMPILING against those packages and
# are never read by an installed SAIGE.so. Not removed and worth knowing why:
# libicudata (30 MB) -- R links ICU for collation and may not start without
# it -- and the ubuntu base (79 MB), whose replacement needs its own testing.
RUN cd .pixi/envs/runtime \
    && rm -rf lib/gcc libexec/gcc x86_64-conda-linux-gnu bin/x86_64-conda-linux-gnu-* \
       lib/libasan* lib/libtsan* lib/libhwasan* lib/liblsan* lib/libubsan* \
       lib/*.a share/doc share/man share/info share/gir-1.0 include \
    && rm -rf lib/R/library/translations conda-meta lib/R/library/*/include \
    && du -sh .

# -------------------------------------------------------------- runtime
FROM ubuntu:24.04 AS runtime
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8 \
    RCPP_PARALLEL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    PATH=/app/.pixi/envs/runtime/bin:/usr/local/bin:/usr/bin:/bin
# Same path as in the build stage: a conda environment is NOT relocatable
# (R_HOME and every activation prefix are baked in at install time).
COPY --from=build /app/.pixi/envs/runtime /app/.pixi/envs/runtime
# pgenlib: harness/build.sh links it with -shared under the name
# libplink2_includes.a, and SAIGE.so records that absolute path (NEEDED), so
# the one file rides along at the same path. ~1 MB.
COPY --from=build /app/.pixi/envs/default/lib/libplink2_includes.a /app/.pixi/envs/default/lib/
COPY --from=build /app/extdata/createSparseGRM.R /app/extdata/step1_fitNULLGLMM.R \
                  /app/extdata/step2_SPAtests.R /app/extdata/extractNglmm.R /usr/local/bin/
# The scripts' shebang runs them through `pixi run` against the repo
# manifest, which serves the dev checkout; there is no pixi here, so it
# becomes the environment's own Rscript (first on PATH above).
RUN sed -i '1s|^#!.*$|#!/usr/bin/env Rscript|' /usr/local/bin/*.R \
    && chmod a+x /usr/local/bin/*.R \
    && createSparseGRM.R --help > /dev/null \
    && step1_fitNULLGLMM.R --help > /dev/null \
    && step2_SPAtests.R --help > /dev/null \
    && extractNglmm.R --help > /dev/null
WORKDIR /data

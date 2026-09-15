# Docker image (strip branch)

A minimal runtime image: R, the packages SAIGE imports, the installed
package, and the four CLI scripts on `PATH`. No compilers, no devtools or
testthat, no plink2, no example data, no goldens, no tests. Mount your own
data and give `--outputPrefix` / `--SAIGEOutputFile` a path on the mount.

## Releases

Pushing a tag `vX.Y.Z` runs `.github/workflows/release.yml`, which refuses
unless `X.Y.Z` equals `DESCRIPTION`'s `Version`, builds `docker/Dockerfile`,
smoke-tests the image (the four `--help`s and one step-2 run on the bundled
100-marker fixture), and pushes

    ghcr.io/astheeggeggs/saige-slim:X.Y.Z
    ghcr.io/astheeggeggs/saige-slim:latest

```sh
git tag v1.5.2 && git push slim v1.5.2
```

## Using it

```sh
docker run --rm -v /path/to/data:/data ghcr.io/astheeggeggs/saige-slim:1.5.2 \
  step1_fitNULLGLMM.R --bedFile=/data/geno.bed --bimFile=/data/geno.bim \
    --famFile=/data/geno.fam --phenoFile=/data/pheno.txt --phenoCol=y \
    --outputPrefix=/data/model
```

The working directory inside the image is `/data`, so relative paths resolve
there when you mount to `/data`. Thread pins (`RCPP_PARALLEL_NUM_THREADS`,
`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`) default to 1 as in
`harness/env.sh`; `--nThreads` overrides the SAIGE ones as usual.

## Building locally

```sh
make docker            # = docker build -f docker/Dockerfile -t saige-slim:dev .
```

Two stages. `build` installs the `default` pixi environment, compiles
pgenlib and `R CMD INSTALL`s the package into the `runtime` environment's
R library (`harness/build.sh` with `SAIGE_R_LIBRARY`). `runtime` starts
from a plain `ubuntu:24.04` and copies only that environment tree and the
scripts. The `runtime` pixi environment is defined in `pixi.toml`
(`[feature.runtime]`, solve group shared with `default` so both resolve
the same R). `.dockerignore` keeps `.pixi/`, `harness/`, `tests/`, `docs/`
and the example data out of the build context.

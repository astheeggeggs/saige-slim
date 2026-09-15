#!/usr/bin/env bash
# Build the documentation: regenerate man/ from the roxygen, then the
# pkgdown site into docs/site/ (gitignored). The tools (roxygen2, pkgdown,
# knitr, rmarkdown, jsonlite, rprojroot) are pinned in the default pixi
# environment:  pixi run bash harness/docs.sh
#
# The article vignettes/articles/cli_options.Rmd calls
# harness/tools/optiondrift.py --json and formals() of the INSTALLED package,
# so build and install first (harness/build.sh). pkgdown articles are not R
# CMD check vignettes (.Rbuildignore excludes vignettes/articles).
#
# Only the rd roclet runs: NAMESPACE is hand-written (exportPattern), and the
# namespace roclet would refuse to touch it and warn on every build.
# load_installed rather than pkgload: pkgload would recompile src/ in-tree,
# which is slow and leaves objects harness/build.sh's incremental guard
# refuses.
#
# CI then runs `git diff --exit-code man/`: a roxygen edit that was not
# followed by this script fails the push. Locally, check `git status man/`.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
# shellcheck disable=SC1091
source harness/env.sh

# The two hand-written pages live in docs/ (many in-repo references point
# there); pkgdown 2.2 reads only .Rmd/.qmd articles from vignettes/articles,
# so each is copied in under a YAML title taken from its first heading (the
# heading itself is dropped so the page does not carry it twice). Neither
# file has a `{r}` chunk, so knitr evaluates nothing. The copies are
# gitignored.
as_article() {  # <source .md> <destination .Rmd>
  local title
  title="$(sed -n 's/^# //p' "$1" | head -1 | sed 's/"/\\"/g')"
  { printf -- '---\ntitle: "%s"\n---\n\n' "$title"; awk 'f || !/^# / { print; next } { f = 1 }' "$1"; } > "$2"
}
as_article docs/cli-options-and-group-file.md vignettes/articles/cli_status.Rmd
as_article docs/flag-surprises.md vignettes/articles/flag_surprises.Rmd

Rscript -e 'roxygen2::roxygenise(".", roclets = "rd", load_code = roxygen2::load_installed)'
Rscript -e 'pkgdown::build_site(".", preview = FALSE, new_process = FALSE, install = FALSE)'
echo "docs: site in docs/site/ ; man/ status:"
git status --short man/

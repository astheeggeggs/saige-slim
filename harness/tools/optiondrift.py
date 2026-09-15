#!/usr/bin/env python3
"""optiondrift.py -- the CLI wrappers' defaults against the exported formals.

  optiondrift.py            report every CLI default that differs from the
                            formal it is forwarded to, minus the allowed list;
                            exit 1 on any unlisted divergence or stale entry
  optiondrift.py selftest   the controls: a planted divergence is caught, a
                            planted agreement is not, and every allowed entry
                            still describes a live divergence
  optiondrift.py --all      also print the agreeing pairs and the unforwarded
                            formals (informational)

Why (LEDGER #151, #115): the four `extdata/*.R` wrappers declare a default per
`make_option()`, the exported function declares another in its formals, and
nothing owned the value -- `--X_PARregion` had THREE (CLI `""`, help text,
formal), `--maxMAF_in_groupTest` gives a CLI user three cells and an API
caller one, `--isCovariateOffset` is TRUE by one route and FALSE by the
other under identical help text. A hand sweep found four; #115's STATUS
records why a hand sweep is not enough. This tool compares VALUES, not
strings: `1L` == `1`, `1e-05` == `1e-5`, a comma string == `c(...)`.

Method, pure Python (no R, no installed package, seconds):
  1. paren-match every `make_option(` block of the wrapper -> (flag, default)
  2. paren-match `<fn> <- function(` in the R source -> (formal, default);
     the step-1 wrapper calls `fitNULLGLMM_multi`, whose `...` forwards to
     `fitNULLGLMM`, so both signatures are read and the multi's own win
  3. paren-match the wrapper's call to <fn> -> `formal = opt$flag`,
     `formal = opt[["flag"]]`, or `formal = <local>` where `<local> <- ...`
     mentions exactly one `opt$flag` (the comma-splitting locals)
  4. parse each default as an R literal (TRUE/FALSE/NULL/NA, numbers, strings,
     `c(...)`); a CLI comma string is split and each piece parsed the same way
     when the formal is a vector or the pieces are all numeric
  5. compare canonical values; a default that is not a literal (a named
     constant such as `flexrvDefaultScoreSpec`) is reported as unparsed and
     never counts as a divergence

An allowed divergence carries the ledger row that records it and BOTH values,
so the entry goes stale -- and the tool fails -- the moment either side moves.
Resolving a divergence means deleting its entry in the same commit.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

# (wrapper script, [(R source, function name)], call-site function name)
WRAPPERS = [
    ("extdata/step2_SPAtests.R",
     [("R/SAIGE_Test_main.R", "SPAGMMATtest")], "SPAGMMATtest"),
    ("extdata/step1_fitNULLGLMM.R",
     [("R/SAIGE_fitGLMM_multi.R", "fitNULLGLMM_multi"),
      ("R/SAIGE_fitGLMM_fast.R", "fitNULLGLMM")], "fitNULLGLMM_multi"),
    ("extdata/createSparseGRM.R",
     [("R/SAIGE_createSparseGRM.R", "createSparseGRM")], "createSparseGRM"),
    ("extdata/extractNglmm.R",
     [("R/SAIGE_extractNeff.R", "getNeff")], "getNeff"),
]

# (wrapper, flag) -> (ledger row, CLI value, formal value). The values are the
# canonical forms this tool prints; an entry whose values no longer match what
# the sources say is STALE and fails the run, because it then allows a
# divergence nobody has looked at.
ALLOWED = {
    # The four the 2026-09-14 hand sweep found (the row's WHERE).
    ("extdata/step2_SPAtests.R", "X_PARregion"):
        (151, "", "60001-2699520,154931044-155270560"),
    ("extdata/step2_SPAtests.R", "maxMAF_in_groupTest"):
        (151, [0.0001, 0.001, 0.01], [0.01]),
    ("extdata/step1_fitNULLGLMM.R", "isCovariateOffset"):
        (151, True, False),
    # The eight this tool found on its first run (2026-09-15, row STATUS).
    ("extdata/step1_fitNULLGLMM.R", "sampleIDColinphenoFile"):
        (151, "IID", ""),
    ("extdata/step1_fitNULLGLMM.R", "outputPrefix"):
        (151, "~/", ""),
    ("extdata/createSparseGRM.R", "nThreads"):
        (151, 16, 1),
    ("extdata/createSparseGRM.R", "outputPrefix"):
        (151, "~/", ""),
    ("extdata/extractNglmm.R", "sampleIDColinphenoFile"):
        (151, "IID", ""),
    ("extdata/extractNglmm.R", "outputPrefix"):
        (151, "~/", ""),
}


# ---------------------------------------------------------------- R lexing

def strip_comments(text):
    """Drop `# ...` to end of line outside string literals."""
    out = []
    i, n = 0, len(text)
    quote = None
    while i < n:
        ch = text[i]
        if quote:
            out.append(ch)
            if ch == "\\" and i + 1 < n:
                out.append(text[i + 1])
                i += 2
                continue
            if ch == quote:
                quote = None
        elif ch in "\"'":
            quote = ch
            out.append(ch)
        elif ch == "#":
            while i < n and text[i] != "\n":
                i += 1
            continue
        else:
            out.append(ch)
        i += 1
    return "".join(out)


def match_paren(text, open_idx):
    """Index of the `)` closing the `(` at open_idx, string-aware."""
    depth = 0
    quote = None
    i = open_idx
    while i < len(text):
        ch = text[i]
        if quote:
            if ch == "\\":
                i += 2
                continue
            if ch == quote:
                quote = None
        elif ch in "\"'":
            quote = ch
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError("unbalanced parentheses")


def split_args(body):
    """Split a call body on top-level commas, string- and paren-aware."""
    parts, depth, quote, cur = [], 0, None, []
    i = 0
    while i < len(body):
        ch = body[i]
        if quote:
            cur.append(ch)
            if ch == "\\" and i + 1 < len(body):
                cur.append(body[i + 1])
                i += 2
                continue
            if ch == quote:
                quote = None
        elif ch in "\"'":
            quote = ch
            cur.append(ch)
        elif ch in "([{":
            depth += 1
            cur.append(ch)
        elif ch in ")]}":
            depth -= 1
            cur.append(ch)
        elif ch == "," and depth == 0:
            parts.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
        i += 1
    tail = "".join(cur).strip()
    if tail:
        parts.append(tail)
    return parts


def balanced_rhs(text, start):
    """The expression starting at `start`, up to the first newline at paren
    depth 0 (a `convertToNumeric(\n  strsplit(...)\n)` spans lines)."""
    depth, quote, i = 0, None, start
    while i < len(text):
        ch = text[i]
        if quote:
            if ch == "\\":
                i += 2
                continue
            if ch == quote:
                quote = None
        elif ch in "\"'":
            quote = ch
        elif ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == "\n" and depth == 0:
            break
        i += 1
    return text[start:i]


def kv(arg):
    """`name = expr` -> (name, expr); bare expr -> (None, expr)."""
    m = re.match(r"^\s*([A-Za-z._][A-Za-z0-9._]*)\s*=(?!=)\s*(.*)$", arg, re.S)
    if m:
        return m.group(1), m.group(2).strip()
    return None, arg.strip()


class Unparsed:
    def __init__(self, text):
        self.text = text

    def __repr__(self):
        return f"<unparsed {self.text!r}>"


NA = "NA"  # sentinel distinct from any string a default could hold


def parse_literal(expr):
    """R literal -> Python value; anything else -> Unparsed."""
    e = expr.strip()
    if e in ("TRUE", "T"):
        return True
    if e in ("FALSE", "F"):
        return False
    if e == "NULL":
        return None
    if e in ("NA", "NA_character_", "NA_real_", "NA_integer_"):
        return NA
    if len(e) >= 2 and e[0] == e[-1] and e[0] in "\"'":
        return e[1:-1]
    m = re.match(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?L?$", e)
    if m:
        return float(e.rstrip("L"))
    m = re.match(r"^c\((.*)\)$", e, re.S)
    if m:
        vals = [parse_literal(p) for p in split_args(m.group(1))]
        if any(isinstance(v, Unparsed) for v in vals):
            return Unparsed(e)
        return vals
    return Unparsed(e)


def canonical(v):
    """Scalars and one-element vectors compare equal; numbers are floats."""
    if isinstance(v, list):
        if len(v) == 1:
            return canonical(v[0])
        return [canonical(x) for x in v]
    return v


def equivalent(cv, fv):
    """Value equality plus the two coercions the R side performs.
    [ASSUMPTION, checked on the current readers: R/PhenoPrep.R reads
    covarColList by length(), R/SAIGE_fitGLMM_multi.R:219 folds a NULL
    outputPrefix_varRatio to "", so NULL == "" == character(0); PhenoPrep
    compares FemaleCode/MaleCode with `==`, and the step-1 wrapper
    as.logical()s --invNormalize, so a string that parses to the formal's
    number or logical is the same value by that route.]"""
    if cv == fv:
        return True
    empty = (None, "", [])
    if cv in empty and fv in empty:
        return True
    if isinstance(cv, str) and isinstance(fv, (bool, float)):
        return parse_literal(cv) == fv
    if isinstance(fv, str) and isinstance(cv, (bool, float)):
        return parse_literal(fv) == cv
    return False


def cli_value(default, formal_value):
    """A CLI default, split on commas when the formal is a vector or every
    piece is numeric (the wrappers' `convertToNumeric(strsplit(...))`)."""
    v = parse_literal(default)
    if isinstance(v, str) and "," in v or (
            isinstance(v, str) and isinstance(formal_value, list)):
        pieces = [p.strip() for p in v.split(",")] if v != "" else []
        parsed = [parse_literal(p) for p in pieces]
        if pieces and all(isinstance(p, float) for p in parsed):
            return parsed
        if isinstance(formal_value, list):
            return pieces
    return v


# ------------------------------------------------------------ extraction

def make_options(text):
    """Every make_option( block -> {flag: default_expr}; a block whose
    `default` is absent gets None, as optparse gives it."""
    text = strip_comments(text)
    out = {}
    for m in re.finditer(r"\bmake_option\s*\(", text):
        close = match_paren(text, m.end() - 1)
        args = split_args(text[m.end():close])
        flag, default, seen_default = None, "NULL", False
        for i, a in enumerate(args):
            k, v = kv(a)
            if k is None and i == 0 or k == "opt_str":
                names = [parse_literal(p) for p in split_args(
                    v[2:-1] if v.startswith("c(") else v)]
                names = [n for n in names if isinstance(n, str)]
                flag = next((n[2:] for n in names if n.startswith("--")), None)
            elif k == "default":
                default, seen_default = v, True
        if flag is None:
            raise ValueError(f"make_option without a --flag: {args[0]!r}")
        out[flag] = default if seen_default else "NULL"
    return out


def formals(text, fn):
    """`fn <- function(...)` -> ordered {name: default_expr or None}."""
    text = strip_comments(text)
    m = re.search(r"(?m)^" + re.escape(fn) + r"\s*<-\s*function\s*\(", text)
    if not m:
        raise ValueError(f"{fn} <- function( not found")
    close = match_paren(text, m.end() - 1)
    out = {}
    for a in split_args(text[m.end():close]):
        k, v = kv(a)
        if k is None:
            out[v] = None  # bare formal (phenoCol, ...)
        else:
            out[k] = v
    return out


def call_map(text, fn):
    """The wrapper's `fn(` call -> {formal: flag} via opt$flag / opt[["flag"]]
    / a local whose assignment mentions exactly one opt$flag."""
    text = strip_comments(text)
    m = re.search(r"(?m)^\s*" + re.escape(fn) + r"\s*\(", text)
    if not m:
        raise ValueError(f"call to {fn}( not found")
    close = match_paren(text, m.end() - 1)
    args = split_args(text[m.end():close])
    # do.call(fn, c(fixed_params, list(...))) is step 2's shape
    if fn == "do.call":
        raise ValueError("do.call form not supported")
    out = {}
    opt_ref = re.compile(r"opt\$([A-Za-z._][A-Za-z0-9._]*)|opt\[\[\"([^\"]+)\"\]\]")
    for a in args:
        k, v = kv(a)
        if k is None:
            continue
        mm = opt_ref.fullmatch(v.strip())
        if mm:
            out[k] = mm.group(1) or mm.group(2)
            continue
        if re.fullmatch(r"[A-Za-z._][A-Za-z0-9._]*", v.strip()):
            local = v.strip()
            am = re.search(r"(?m)^" + re.escape(local) + r"\s*<-\s*", text)
            if am:
                rhs = balanced_rhs(text, am.end())
                refs = {g1 or g2 for g1, g2 in opt_ref.findall(rhs)}
                if len(refs) == 1:
                    out[k] = refs.pop()
    return out


def step2_call_map(text):
    """Step 2 builds `fixed_params <- list(...)` then do.call()s it."""
    text = strip_comments(text)
    m = re.search(r"(?m)^fixed_params\s*<-\s*list\s*\(", text)
    close = match_paren(text, m.end() - 1)
    fake = "SPAGMMATtest(" + text[m.end():close] + ",\n" + \
        "SAIGEOutputFile = opt$SAIGEOutputFile, idstoIncludeFile = opt$idstoIncludeFile, " \
        "nThreads = opt$nThreads)"
    # keep the local assignments visible for the `<local> <- ...` lookup
    return call_map(text[:m.start()] + "\n" + fake, "SPAGMMATtest")


# ------------------------------------------------------------- comparison

def compare(wrapper, cli, forms, cmap):
    """-> (divergences, agreements, unparsed, unforwarded). A divergence is
    (flag, formal, cli_value, formal_value)."""
    div, agree, unparsed = [], [], []
    for formal, flag in cmap.items():
        if formal not in forms or flag not in cli:
            continue
        fexpr = forms[formal]
        if fexpr is None:
            continue  # bare formal: no default to disagree with
        fv = parse_literal(fexpr)
        if isinstance(fv, Unparsed):
            unparsed.append((flag, formal, fexpr))
            continue
        cv = cli_value(cli[flag], fv)
        if isinstance(cv, Unparsed):
            unparsed.append((flag, formal, cli[flag]))
            continue
        cv, fv = canonical(cv), canonical(fv)
        (div if not equivalent(cv, fv) else agree).append((flag, formal, cv, fv))
    unforwarded = sorted(set(forms) - set(cmap))
    return div, agree, unparsed, unforwarded


def load(wrapper, sources, callfn):
    wtext = open(os.path.join(REPO, wrapper)).read()
    cli = make_options(wtext)
    forms = {}
    for src, fn in reversed(sources):  # first-listed wins on overlap
        forms.update(formals(open(os.path.join(REPO, src)).read(), fn))
    cmap = step2_call_map(wtext) if callfn == "SPAGMMATtest" else call_map(wtext, callfn)
    return cli, forms, cmap


def check(show_all=False):
    rc = 0
    seen_allowed = set()
    for wrapper, sources, callfn in WRAPPERS:
        cli, forms, cmap = load(wrapper, sources, callfn)
        div, agree, unparsed, unforwarded = compare(wrapper, cli, forms, cmap)
        for flag, formal, cv, fv in div:
            key = (wrapper, flag)
            if key in ALLOWED:
                row, acv, afv = ALLOWED[key]
                seen_allowed.add(key)
                if canonical(acv) != cv or canonical(afv) != fv:
                    print(f"STALE ALLOWED (#{row}) {wrapper} --{flag}: entry says "
                          f"CLI {acv!r} / formal {afv!r}, sources say {cv!r} / {fv!r}")
                    rc = 1
                elif show_all:
                    print(f"allowed (#{row}) {wrapper} --{flag}: CLI {cv!r} / formal {fv!r}")
                continue
            print(f"DRIFT {wrapper} --{flag} -> {formal}: CLI {cv!r} / formal {fv!r}")
            rc = 1
        if show_all:
            for flag, formal, cv, fv in agree:
                print(f"ok    {wrapper} --{flag} -> {formal}: {cv!r}")
            for flag, formal, expr in unparsed:
                print(f"skip  {wrapper} --{flag} -> {formal}: not a literal: {expr}")
            for formal in unforwarded:
                print(f"info  {wrapper}: formal {formal} not forwarded by the CLI")
        print(f"{wrapper}: {len(cmap)} forwarded, {len(agree)} agree, "
              f"{len(div)} diverge, {len(unparsed)} unparsed, {len(unforwarded)} unforwarded")
    for key in set(ALLOWED) - seen_allowed:
        print(f"STALE ALLOWED (#{ALLOWED[key][0]}) {key[0]} --{key[1]}: no longer diverges; "
              f"delete the entry")
        rc = 1
    return rc


# --------------------------------------------------------------- controls

def selftest():
    # Positive control: a planted divergence in a synthetic wrapper is caught.
    wrapper = '''
option_list <- list(
  make_option("--a", type = "numeric", default = 1L, help = "x"),
  make_option("--b", type = "character", default = "0.1,0.2", help = "y # not a comment"),
  make_option("--c", type = "logical", default = TRUE, help = "z"),
  make_option("--d", type = "numeric", default = 1e-05, help = "w")
)
b_local <- as.numeric(strsplit(opt$b, ",")[[1]])
fn(
  a = opt$a,
  b = b_local,
  c = opt[["c"]],
  d = opt$d
)
'''
    source = '''
fn <- function(a = 1, # comment with ) paren
               b = c(0.1, 0.2),
               c = FALSE,
               d = 1e-5) {
}
'''
    div, agree, unparsed, _ = compare(
        "synthetic", make_options(wrapper), formals(source, "fn"), call_map(wrapper, "fn"))
    got = {d[0] for d in div}
    ok = {a[0] for a in agree}
    if got != {"c"} or ok != {"a", "b", "d"} or unparsed:
        print(f"optiondrift selftest FAILED: diverge={got} agree={ok} unparsed={unparsed}")
        return 1
    # Negative control on the real sources: every ALLOWED entry names a
    # divergence that is still live (check() reports it stale otherwise).
    rc = check(show_all=False)
    if rc:
        print("optiondrift selftest FAILED: the real sources do not match ALLOWED")
        return 1
    print(f"optiondrift selftest OK: planted divergence caught, planted agreements "
          f"(1L/1, comma-string/c(), 1e-05/1e-5) pass, {len(ALLOWED)} allowed entries live")
    return 0


def main(argv):
    if argv and argv[0] == "selftest":
        return selftest()
    return check(show_all="--all" in argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

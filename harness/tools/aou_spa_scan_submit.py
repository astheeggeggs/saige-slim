# Notebook cell: launch harness/tools/aou_spa_scan.sh over one delivered output
# directory as a dsub job (uses the notebook's aou_dsub, BUCKET, RUNS conventions).
# One job per outdir; the FlexRV arm of one ancestry streams ~25 GB and takes
# under an hour on n2-standard-2 -- nothing is localized, nothing egresses.
import shlex, subprocess, os, re
from datetime import datetime

def submit_scan(outdir, pattern="*", tail_p="1e-5", image="google/cloud-sdk:slim",
                machine_type="n2-standard-2", dry=False):
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    scan_dir = f"{BUCKET}/SAIGE_SLIM_TESTING/{outdir}"
    out_dir  = f"{BUCKET}/SAIGE_SLIM_TESTING/spa_scan/{outdir}/{stamp}"
    name = "spascan-" + re.sub(r"[._]", "-", outdir).lower()
    args = ["aou_dsub", "--name", name, "--image", image,
            "--env", f"SCAN_DIR={scan_dir}", "--env", f"TAIL_P={tail_p}",
            "--env", f"PATTERN={pattern}",
            "--output-recursive", f"OUT_DIR={out_dir}",
            "--machine-type", machine_type, "--disk-size", "20",
            "--script", "aou_spa_scan.sh"]
    if dry: args.insert(1, "--dry-run")
    script = "source /home/jupyter/aou_dsub.bash\n" + " ".join(shlex.quote(a) for a in args)
    r = subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=os.environ)
    print(r.stdout[-1500:] or r.stderr[-1500:])
    print("results will land in", out_dir)
    return out_dir

# e.g. one per delivered directory; PATTERN narrows to an arm ("fxv4.*", "skv4.*")
# submit_scan("<outdir>", pattern="fxv4.*")

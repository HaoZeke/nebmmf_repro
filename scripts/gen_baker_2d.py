#!/usr/bin/env python3
"""Generate all 48 baker_2d landscape plots in parallel with xyzrender."""
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

PAPER = "/home/goswami/Git/Github/TeX/roneb_tex/v2/FrontiersChemistry/imgs"
NEB = "results/03_neb"
PYTHON = sys.executable
BASE_CWD = os.getcwd()
CALLS_RE = re.compile(r"(\d+) total_force_calls")

REACTIONS = {
    "01_hcn": "HCN -> HNC",
    "02_hcch": "HCCH -> CCH2",
    "03_h2co": "H2CO -> H2 + CO",
    "04_ch3o": "CH3O -> CH2OH",
    "05_cyclopropyl": "Cyclopropyl ring opening",
    "06_bicyclobutane": "Bicyclobutane -> butadiene",
    "08_formyloxyethyl": "Formyloxyethyl 1,2-migration",
    "09_parentdielsalder": "Diels-Alder",
    "10_tetrazine": "s-Tetrazine -> 2HCN + N2",
    "11_trans_butadiene": "trans- -> cis-butadiene",
    "12_ethane_h2_abstraction": "CH3CH3 -> CH2CH2 + H2",
    "13_hf_abstraction": "CH3CH2F -> CH2CH2 + HF",
    "14_vinyl_alcohol": "Acetaldehyde keto-enol",
    "15_hocl": "HCOCl -> HCl + CO",
    "16_h2po4_anion": "H2O + PO3- -> H2PO4-",
    "17_claisen": "Claisen rearrangement",
    "18_silylene_insertion": "SiH2 + CH3CH3",
    "19_hnccs": "HNCCS -> HNC + CS",
    "20_hconh3_cation": "HCONH3+ -> NH4+ + CO",
    "21_acrolein_rot": "Acrolein rotational TS",
    "22_hconhoh": "HCONHOH -> HCOHNHO",
    "23_hcn_h2": "HNC + H2 -> H2CNH",
    "24_h2cnh": "H2CNH -> HCNH2",
    "25_hcnh2": "HCNH2 -> HCN + H2",
}

METHOD_LABELS = {"cineb": "CI-NEB", "mmf": "OCI-NEB"}


def run_plot(args):
    sys_name, method = args
    other = "mmf" if method == "cineb" else "cineb"
    base = f"{BASE_CWD}/{NEB}/{sys_name}/{method}"
    con = f"{base}/neb.con"
    dat_pattern = f"{base}/neb_*.dat"
    path_pattern = f"{base}/neb_path*.con"
    out = f"{PAPER}/suppl/baker_2d/{sys_name}/{method}_landscape.png"

    # Get force calls for title
    results_dat = Path(f"{base}/results.dat")
    calls = "?"
    if results_dat.exists():
        m = CALLS_RE.search(results_dat.read_text())
        if m:
            calls = m.group(1)

    sys_num = sys_name.split("_")[0]
    reaction = REACTIONS.get(sys_name, sys_name)
    method_label = METHOD_LABELS.get(method, method)
    title = f"{sys_num}: {reaction} ({method_label}, {calls} calls)"

    cmd = [
        PYTHON, "-m", "rgpycrumbs.cli", "--dev", "eon", "plt-neb",
        "--con-file", con,
        "--output-file", out,
        "--plot-type", "landscape",
        "--rc-mode", "rmsd",
        "--surface-type", "grad_imq",
        "--show-pts",
        "--landscape-path", "all",
        "--input-dat-pattern", dat_pattern,
        "--input-path-pattern", path_pattern,
        "--ira-kmax", "14",
        "--plot-structures", "crit_points",
        "--dpi", "300",
        "--figsize", "5.37", "5.37",
        "--fontsize-base", "12",
        "--show-legend",
        "--rotation", "auto",
        "--strip-renderer", "xyzrender",
        "--zoom-ratio", "0.3",
        "--strip-dividers",
        "--title", title,
    ]

    # Per-system GP overrides for systems where grad_imq fails
    SURFACE_OVERRIDES = {("06_bicyclobutane", "cineb"): "grad_matern"}
    surface_override = SURFACE_OVERRIDES.get((sys_name, method))
    if surface_override:
        # Replace the surface type in the cmd
        st_idx = cmd.index("--surface-type") + 1
        cmd[st_idx] = surface_override

    # Overlay other method's saddle point
    sp_con = f"{BASE_CWD}/{NEB}/{sys_name}/{other}/sp.con"
    if Path(sp_con).exists():
        cmd.extend(["--additional-con", sp_con, other.upper()])

    # Run in a per-worker temp directory so .neb_landscape.parquet
    # caches don't collide between workers
    with tempfile.TemporaryDirectory(prefix=f"baker2d_{sys_name}_{method}_") as tmpdir:
        r = subprocess.run(cmd, capture_output=True, timeout=180, cwd=tmpdir)

    if r.returncode != 0:
        err = r.stderr.decode()[-300:] if r.stderr else "unknown"
        last = err.strip().splitlines()[-1] if err.strip() else "?"
        return f"{sys_name}/{method}: FAIL ({last})"
    return f"{sys_name}/{method}: OK"



def run_1d_profile(args):
    sys_name, method = args
    base = f"{BASE_CWD}/{NEB}/{sys_name}/{method}"
    con = f"{base}/neb.con"
    dat_pattern = f"{base}/neb_*.dat"
    out = f"{PAPER}/suppl/baker_2d/{sys_name}/{method}_1D_path.png"

    results_dat = Path(f"{base}/results.dat")
    calls = "?"
    if results_dat.exists():
        m = CALLS_RE.search(results_dat.read_text())
        if m:
            calls = m.group(1)

    sys_num = sys_name.split("_")[0]
    reaction = REACTIONS.get(sys_name, sys_name)
    method_label = METHOD_LABELS.get(method, method)
    title = f"{sys_num}: {reaction} ({method_label}, {calls} calls)"

    cmd = [
        PYTHON, "-m", "rgpycrumbs.cli", "--dev", "eon", "plt-neb",
        "--con-file", con,
        "--output-file", out,
        "--plot-type", "profile",
        "--rc-mode", "path",
        "--input-dat-pattern", dat_pattern,
        "--plot-structures", "crit_points",
        "--dpi", "300",
        "--figsize", "5.37", "3.37",
        "--fontsize-base", "12",
        "--rotation", "auto",
        "--strip-renderer", "xyzrender",
        "--zoom-ratio", "0.1",
        "--strip-dividers",
        "--no-legend",
        "--title", title,
    ]

    r = subprocess.run(cmd, capture_output=True, timeout=180)
    if r.returncode != 0:
        err = r.stderr.decode()[-300:] if r.stderr else "unknown"
        last = err.strip().splitlines()[-1] if err.strip() else "?"
        return f"{sys_name}/{method} 1D: FAIL ({last})"
    return f"{sys_name}/{method} 1D: OK"


if __name__ == "__main__":
    systems = sorted(d.name for d in Path(NEB).iterdir() if d.is_dir())
    tasks = [(s, m) for s in systems for m in ["cineb", "mmf"]]

    for s in systems:
        Path(f"{PAPER}/suppl/baker_2d/{s}").mkdir(parents=True, exist_ok=True)

    # Delete any stale root cache
    stale = Path(".neb_landscape.parquet")
    if stale.exists():
        stale.unlink()

    print(f"Generating {len(tasks)} plots on 16 cores...")
    ok = fail = 0
    with ProcessPoolExecutor(max_workers=16) as ex:
        futures = {ex.submit(run_plot, t): t for t in tasks}
        for f in as_completed(futures):
            result = f.result()
            print(result, flush=True)
            if "OK" in result:
                ok += 1
            else:
                fail += 1
    # 1D profiles
    print(f"\nGenerating {len(tasks)} 1D profiles...")
    with ProcessPoolExecutor(max_workers=16) as ex:
        futures_1d = {ex.submit(run_1d_profile, t): t for t in tasks}
        for ft in as_completed(futures_1d):
            result = ft.result()
            print(result, flush=True)
            if "OK" in result:
                ok += 1
            else:
                fail += 1
    print(f"\nTotal: {ok} OK, {fail} FAIL")

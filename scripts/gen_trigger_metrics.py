#!/usr/bin/env python3
"""Extract OCI-NEB trigger/backoff metrics from client_quill.log files.

Parses MMF trigger and backoff events from eOn NEB logs and outputs
an org-mode table suitable for the supplementary material.

Usage:
    python scripts/gen_trigger_metrics.py [NEB_DIR]

NEB_DIR defaults to eonRuns/results/03_neb.
"""

import re
import sys
from pathlib import Path

BAKER_SYSTEMS = [
    "01_hcn", "02_hcch", "03_h2co", "04_ch3o", "05_cyclopropyl",
    "06_bicyclobutane", "08_formyloxyethyl", "09_parentdielsalder",
    "10_tetrazine", "11_trans_butadiene", "12_ethane_h2_abstraction",
    "13_hf_abstraction", "14_vinyl_alcohol", "15_hocl", "16_h2po4_anion",
    "17_claisen", "18_silylene_insertion", "19_hnccs", "20_hconh3_cation",
    "21_acrolein_rot", "22_hconhoh", "23_hcn_h2", "24_h2cnh", "25_hcnh2",
]

TRIGGER_RE = re.compile(r"Triggering MMF")
BACKOFF_RE = re.compile(r"(?:MMF |)[Bb]ackoff")
CONVERGED_MMF_RE = re.compile(r"(?:NEB )?converged after MMF|Converged.*MMF")
CONVERGED_NEB_RE = re.compile(r"NEB converged|converged_force")


def parse_log(log_path: Path) -> tuple[int, int, str]:
    """Return (n_triggers, n_backoffs, final_state) from a client_quill.log."""
    if not log_path.exists():
        return 0, 0, "Missing"

    txt = log_path.read_text(errors="ignore")
    n_triggers = len(TRIGGER_RE.findall(txt))
    n_backoffs = len(BACKOFF_RE.findall(txt))

    if CONVERGED_MMF_RE.search(txt):
        state = "Converged (MMF)"
    elif CONVERGED_NEB_RE.search(txt):
        state = "Converged (NEB)"
    else:
        state = "Unknown"

    return n_triggers, n_backoffs, state


def main():
    script_dir = Path(__file__).resolve().parent
    base_dir = script_dir.parent / "eonRuns"
    if not base_dir.is_dir():
        base_dir = Path.cwd()

    neb_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else base_dir / "results" / "03_neb"

    print("| system                   | n_triggers | n_backoffs | final_state     |")
    print("|--------------------------+------------+------------+-----------------|")

    for sys_name in BAKER_SYSTEMS:
        log = neb_dir / sys_name / "mmf" / "client_quill.log"
        n_trig, n_back, state = parse_log(log)
        print(f"| {sys_name:<24s} | {n_trig:>10d} | {n_back:>10d} | {state:<15s} |")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "rich",
#   "click",
#   "polars",
#   "ase",
#   "numpy",
#   "rgpycrumbs @ git+https://github.com/HaoZeke/rgpycrumbs@eonMLFlow",
# ]
# ///

"""
MMF/Dimer Analysis Tool.

Analyzes single-ended transition state searches (e.g., from '04_mmf') by:
1. Parsing 'results.dat' for computational cost (time, calls) and energetics.
2. Parsing 'saddle.con' for the final transition state geometry.
3. Parsing 'pos.con' (if available) to calculate the "Dimer Correction Distance" (RMSD travelled from seed).
4. Reporting a consolidated table of results.
"""

import logging
import pathlib
import re
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from enum import IntEnum

import click
import numpy as np
import polars as pl
from ase import Atoms
from ase.io import read as ase_read
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress

# --- Import User's Alignment Module ---
try:
    from rgpycrumbs.geom.api.alignment import (
        IRAConfig,
        align_structure_robust,
    )
except ImportError:
    align_structure_robust = None

# --- Configuration ---

logging.basicConfig(
    level="WARNING",
    format="%(message)s",
    handlers=[RichHandler(show_path=False)],
)

# Regex patterns
CALLS_REGEX = re.compile(r"destroyed after (\d+) calls")
# Matches typical EON output for final energy or extrema
ENERGY_REGEX = re.compile(r"energy\s+([\d.-]+)")
TIME_REGEX = re.compile(r"time_seconds\s+([\d.]+)")


class RunStatus(IntEnum):
    """Maps EON termination codes to human-readable statuses."""

    GOOD = 0
    INIT = 1
    BAD_MAX_ITERATIONS = 2
    RUNNING = 3
    MAX_UNCERTAINTY = 4
    UNKNOWN = 999

    @classmethod
    def from_code(cls, code: str) -> "RunStatus":
        try:
            return cls(int(code))
        except (ValueError, TypeError):
            return cls.UNKNOWN


@dataclass
class RunMetrics:
    system: str
    method: str
    time_seconds: float | None = None
    calls: int | None = None
    potential_energy: float | None = None
    max_force: float | None = None

    # Geometry
    saddle_atoms: Atoms | None = None

    # Distance from Initial Guess (pos.con) to Final Saddle (saddle.con)
    dimer_travel_dist: float | None = None

    status: RunStatus = RunStatus.UNKNOWN


REACTION_MAP = {
    "01": "bold('01:')~HCN*'->'*HNC",
    "02": "bold('02:')~HCCH*'->'*CCH[2]",
    "03": "bold('03:')~H[2]*CO*'->'*H[2] + CO",
    "04": "bold('04:')~CH[3]*O*'->'*CH[2]*OH",
    "05": "bold('05:')~Cyclopropyl~ring~opening",
    "06": "bold('06:')~Bicyclo1.1.0*butane*'->'*italic(trans)-butadiene",
    "08": "bold('08:')~Formyloxyethyl~1*','*2-migration",
    "09": "bold('09:')~Parent~Diels-Alder~cycloaddition",
    "10": "bold('10:')~italic(s)-tetrazine*'->'*2*HCN + N[2]",
    "11": "bold('11:')~italic(trans)-butadiene*'->'*italic(cis)-butadiene",
    "12": "bold('12:')~CH[3]*CH[3]*'->'*CH[2]*CH[2] + H[2]",
    "13": "bold('13:')~CH[3]*CH[2]*F*'->'*CH[2]*CH[2] + HF",
    "14": "bold('14:')~Acetaldehyde~keto-enol~tautomerism",
    "15": "bold('15:')~HOCl*'->'*HCl + CO",
    "16": "bold('16:')~H[2]*O + PO[3]^-{}*'->'*H[2]*PO[4]^-{}",
    "17": "bold('17:')~CH[2]*CHCH[2]*CH[2]*CHO~Claisen~rearrangement",
    "18": "bold('18:')~SiH[2] + CH[3]*CH[3]*'->'*SiH[3]*CH[2]*CH[3]",
    "19": "bold('19:')~HNCCS*'->'*HNC + CS",
    "20": "bold('20:')~HCONH[3]^+{}*'->'*NH[4]^+{} + CO",
    "21": "bold('21:')~Acrolein~rotational~TS",
    "22": "bold('22:')~HCONHOH*'->'*HCOHNHO",
    "23": "bold('23:')~HNC + H[2]*'->'*H[2]*CNH",
    "24": "bold('24:')~H[2]*CNH*'->'*HCNH[2]",
    "25": "bold('25:')~HCNH[2]*'->'*HCN + H[2]",
}


def calculate_rmsd(atoms_a: Atoms, atoms_b: Atoms) -> float:
    """Calculates RMSD using robust alignment (IRA/Procrustes)."""
    if not align_structure_robust or not atoms_a or not atoms_b:
        return -1.0

    try:
        mobile = atoms_b.copy()
        ref = atoms_a.copy()
        config = IRAConfig(enabled=True, kmax=1.8)

        # Aligns 'mobile' onto 'ref' in-place
        align_structure_robust(ref, mobile, config)

        pos_a = ref.get_positions()
        pos_b = mobile.get_positions()
        diff = pos_a - pos_b
        rmsd = np.sqrt((diff * diff).sum() / len(ref))
        return rmsd
    except Exception:
        return -1.0


def parse_run(path: pathlib.Path) -> RunMetrics | None:
    """Parses a single run directory."""

    # 1. Deduce System/Method from Path
    # Assumption: .../system_name/method_name/results.dat
    try:
        # parent is method, parent.parent is system
        method_name = path.parent.name
        system_name = path.parent.parent.name
    except IndexError:
        return None

    results_file = path
    run_dir = path.parent

    # 2. Parse results.dat
    time_val, energy_val = None, None
    status = RunStatus.UNKNOWN

    if results_file.exists():
        txt = results_file.read_text(errors="ignore")

        # Time
        t_match = TIME_REGEX.search(txt)
        if t_match:
            time_val = float(t_match.group(1))

        # Energy
        # EON results often list the final potential energy
        e_match = ENERGY_REGEX.search(txt)
        if e_match:
            energy_val = float(e_match.group(1))

        # Termination Status
        first_line = txt.splitlines()[0] if txt else ""
        if "termination_reason" in first_line:
            parts = first_line.split()
            if parts:
                status = RunStatus.from_code(parts[0])

    # 3. Parse Calls from logs
    calls = None
    for log in run_dir.glob("*.log"):
        try:
            l_txt = log.read_text(errors="ignore")
            c_match = CALLS_REGEX.search(l_txt)
            if c_match:
                calls = int(c_match.group(1))
                break
        except Exception:
            pass

    # 4. Load Geometries
    saddle_atoms = None
    start_atoms = None

    # Saddle (Result)
    saddle_path = run_dir / "saddle.con"
    if not saddle_path.exists():
        # Fallback to min.con if saddle.con wasn't renamed
        saddle_path = run_dir / "min.con"

    if saddle_path.exists():
        try:
            saddle_atoms = ase_read(saddle_path)
        except Exception:
            pass

    # Start (Seed)
    pos_path = run_dir / "pos.con"
    if pos_path.exists():
        try:
            start_atoms = ase_read(pos_path)
        except Exception:
            pass

    # 5. Calculate Travel Distance (Dimer Correction)
    dist = None
    if saddle_atoms and start_atoms:
        if len(saddle_atoms) == len(start_atoms):
            dist = calculate_rmsd(start_atoms, saddle_atoms)

    # 6. Extract Max Force if available (often in min.con or results)
    # Simple heuristic: calculate from atoms if calculator not attached
    # Here we skip unless we parse force explicitly.
    # For now, we leave max_force None or assume parsed from log if needed.

    return RunMetrics(
        system=system_name,
        method=method_name,
        time_seconds=time_val,
        calls=calls,
        potential_energy=energy_val,
        saddle_atoms=saddle_atoms,
        dimer_travel_dist=dist,
        status=status,
    )


@click.command()
@click.argument(
    "search_path", type=click.Path(exists=True, path_type=pathlib.Path), default="."
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=pathlib.Path),
    default="neb_mmf_static_results.csv",
    help="Output CSV path.",
)
def main(search_path: pathlib.Path, output: pathlib.Path):
    """
    Crawls SEARCH_PATH for 'results.dat' files, assumes MMF/Dimer structure,
    and tabulates performance metrics and saddle point geometries.
    """
    console = Console()

    # 1. Find runs
    results_files = list(search_path.rglob("results.dat"))
    if not results_files:
        console.print("[red]No results.dat files found![/]")
        return

    console.print(f"[bold green]Found {len(results_files)} runs. Parsing...[/]")

    # 2. Parse in Parallel
    runs = []
    with ThreadPoolExecutor() as exc:
        futures = [exc.submit(parse_run, f) for f in results_files]
        for f in futures:
            res = f.result()
            if res:
                runs.append(res)

    if not runs:
        console.print("[yellow]No valid run data extracted.[/]")
        return

    # 3. Build DataFrame
    data = []
    with Progress(console=console) as progress:
        task = progress.add_task("[cyan]Processing Geometries...", total=len(runs))

        for r in runs:
            # Metadata
            formula = "Unknown"
            mass = 0.0
            natoms = 0

            if r.saddle_atoms:
                formula = r.saddle_atoms.get_chemical_formula()
                mass = sum(r.saddle_atoms.get_masses())
                natoms = len(r.saddle_atoms)

            # Label
            # Try to grab the number prefix from system name (e.g. "01_HCN...")
            rxn_label = r.system
            prefix_match = re.match(r"^(\d+)", r.system)
            if prefix_match:
                key = prefix_match.group(1)
                if key in REACTION_MAP:
                    rxn_label = REACTION_MAP[key]

            data.append(
                {
                    "System": r.system,
                    "Method": r.method,
                    "Reaction": rxn_label,
                    "Formula": formula,
                    "N_Atoms": natoms,
                    "Mass": mass,
                    "Energy": r.potential_energy,
                    "Time_s": r.time_seconds,
                    "Calls": r.calls,
                    "Dimer_Travel_RMSD": r.dimer_travel_dist,
                    "Status": r.status.name,
                }
            )
            progress.advance(task)

    df = pl.DataFrame(data)

    # 4. Save
    # Sort by complexity then system
    df = df.sort(["N_Atoms", "System"])
    df.write_csv(output)

    # 5. Display Table
    console.print("\n[bold]Run Summary[/]")

    # Select columns for display
    display_cols = [
        "System",
        "Formula",
        "Time_s",
        "Calls",
        "Dimer_Travel_RMSD",
        "Status",
    ]

    # Filter to existing columns (just in case)
    display_cols = [c for c in display_cols if c in df.columns]

    pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=console.width)
    console.print(df.select(display_cols))

    console.print(f"\n[bold green]Saved results to {output}[/]")


if __name__ == "__main__":
    main()

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
Combined NEB + MMF Benchmark Tool.

Analyzes the full workflow (NEB Pre-stage + MMF/Dimer Final Stage) by:
1. Crawling for 'results.dat' files.
2. Grouping runs by System and Method.
3. Calculating cumulative metrics (Total Time, Total Force Calls).
4. Analyzing the final saddle point geometry and the distance traveled during the Dimer stage.
"""

import logging
import pathlib
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
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
class SingleStageMetrics:
    """Metrics for one stage (either NEB or MMF)."""

    time_seconds: float = 0.0
    calls: int = 0
    energy: float | None = None
    status: RunStatus = RunStatus.UNKNOWN
    path: pathlib.Path | None = None


@dataclass
class CombinedMetrics:
    """Aggregated metrics for a System/Method pair."""

    system: str
    method: str

    # Per-stage data
    neb: SingleStageMetrics = field(default_factory=SingleStageMetrics)
    mmf: SingleStageMetrics = field(default_factory=SingleStageMetrics)

    # Geometry (Final Saddle from MMF)
    saddle_atoms: Atoms | None = None

    # Analysis
    dimer_travel_rmsd: float | None = (
        None  # Distance from MMF start (NEB peak) to MMF finish
    )

    @property
    def total_time(self) -> float:
        return self.neb.time_seconds + self.mmf.time_seconds

    @property
    def total_calls(self) -> int:
        return self.neb.calls + self.mmf.calls


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
    """Calculates RMSD using robust alignment."""
    if not align_structure_robust or not atoms_a or not atoms_b:
        return -1.0
    try:
        mobile = atoms_b.copy()
        ref = atoms_a.copy()
        config = IRAConfig(enabled=True, kmax=1.8)
        align_structure_robust(ref, mobile, config)
        diff = ref.get_positions() - mobile.get_positions()
        return np.sqrt((diff**2).sum() / len(ref))
    except Exception:
        return -1.0


def parse_stage_metrics(results_path: pathlib.Path) -> SingleStageMetrics:
    """Extracts basic metrics from a results.dat file."""
    metrics = SingleStageMetrics(path=results_path.parent)

    if results_path.exists():
        txt = results_path.read_text(errors="ignore")

        # Time
        t_match = TIME_REGEX.search(txt)
        if t_match:
            metrics.time_seconds = float(t_match.group(1))

        # Energy
        e_match = ENERGY_REGEX.search(txt)
        if e_match:
            metrics.energy = float(e_match.group(1))

        # Status
        first_line = txt.splitlines()[0] if txt else ""
        if "termination_reason" in first_line:
            parts = first_line.split()
            if parts:
                metrics.status = RunStatus.from_code(parts[0])

    # Calls from logs
    for log in results_path.parent.glob("*.log"):
        try:
            l_txt = log.read_text(errors="ignore")
            c_match = CALLS_REGEX.search(l_txt)
            if c_match:
                metrics.calls = int(c_match.group(1))
                break
        except Exception:
            pass

    return metrics


def process_single_file(path: pathlib.Path):
    """Worker function to parse a file and identify its context."""
    # Heuristic: .../{stage}/{system}/{method}/results.dat
    # stages usually: 03_neb (or similar containing 'neb') / 04_mmf (or similar containing 'mmf')

    try:
        parts = path.parts
        # Look specifically for parent folders named roughly 'neb' or 'mmf'
        method_idx = -2
        system_idx = -3
        stage_idx = -4

        stage_dir = parts[stage_idx]
        system_name = parts[system_idx]
        method_name = parts[method_idx]

        stage_type = "unknown"
        if "neb" in stage_dir.lower():
            stage_type = "neb"
        elif "mmf" in stage_dir.lower():
            stage_type = "mmf"

        metrics = parse_stage_metrics(path)
        return system_name, method_name, stage_type, metrics

    except IndexError:
        return None


@click.command()
@click.argument(
    "search_path", type=click.Path(exists=True, path_type=pathlib.Path), default="."
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=pathlib.Path),
    default="combined_static_results.csv",
    help="Output CSV path.",
)
def main(search_path: pathlib.Path, output: pathlib.Path):
    """
    Analyzes NEB+MMF workflows.
    Calculates Total Calls and Total Time across both stages.
    """
    console = Console()

    # 1. Find all results
    results_files = list(search_path.rglob("results.dat"))
    if not results_files:
        console.print("[red]No results.dat files found![/]")
        return

    console.print(f"[bold green]Parsing {len(results_files)} output files...[/]")

    # 2. Parse and Group
    # Map: (system, method) -> CombinedMetrics
    grouped_data = defaultdict(lambda: {"neb": None, "mmf": None})
    systems_seen = set()
    methods_seen = set()

    with ThreadPoolExecutor() as exc:
        futures = [exc.submit(process_single_file, p) for p in results_files]
        for f in futures:
            res = f.result()
            if res:
                sys, meth, stage, metrics = res
                if stage in ["neb", "mmf"]:
                    grouped_data[(sys, meth)][stage] = metrics
                    systems_seen.add(sys)
                    methods_seen.add(meth)

    # 3. Aggregate and Analyze Geometry
    final_results = []

    with Progress(console=console) as progress:
        task = progress.add_task(
            "[cyan]Analyzing Combined Workflows...", total=len(grouped_data)
        )

        for (sys, meth), stages in grouped_data.items():
            combined = CombinedMetrics(system=sys, method=meth)

            # Populate Stages
            if stages["neb"]:
                combined.neb = stages["neb"]
            if stages["mmf"]:
                combined.mmf = stages["mmf"]

            # Geometry Analysis (Focused on MMF stage)
            if combined.mmf.path:
                # Load Saddle
                saddle_p = combined.mmf.path / "saddle.con"
                if not saddle_p.exists():
                    saddle_p = combined.mmf.path / "min.con"

                if saddle_p.exists():
                    try:
                        atoms = ase_read(saddle_p)
                        combined.saddle_atoms = atoms
                    except Exception:
                        pass

                # Load Seed (pos.con) to calc travel
                start_p = combined.mmf.path / "pos.con"
                if start_p.exists() and combined.saddle_atoms:
                    try:
                        start_atoms = ase_read(start_p)
                        combined.dimer_travel_rmsd = calculate_rmsd(
                            start_atoms, combined.saddle_atoms
                        )
                    except Exception:
                        pass

            final_results.append(combined)
            progress.advance(task)

    # 4. Export to CSV
    csv_rows = []
    for r in final_results:
        # Metadata
        natoms = len(r.saddle_atoms) if r.saddle_atoms else 0
        formula = r.saddle_atoms.get_chemical_formula() if r.saddle_atoms else "Unknown"

        # Reaction Label
        rxn_label = r.system
        prefix_match = re.match(r"^(\d+)", r.system)
        if prefix_match and prefix_match.group(1) in REACTION_MAP:
            rxn_label = REACTION_MAP[prefix_match.group(1)]

        csv_rows.append(
            {
                "System": r.system,
                "Method": r.method,
                "Reaction": rxn_label,
                "N_Atoms": natoms,
                "Formula": formula,
                # Totals
                "Total_Time_s": r.total_time,
                "Total_Calls": r.total_calls,
                # MMF Specifics (Final Result)
                "Final_Energy": r.mmf.energy,
                "MMF_Status": r.mmf.status.name,
                "Dimer_Travel_RMSD": r.dimer_travel_rmsd,
                "MMF_Time": r.mmf.time_seconds,
                "MMF_Calls": r.mmf.calls,
                # NEB Specifics
                "NEB_Time": r.neb.time_seconds,
                "NEB_Calls": r.neb.calls,
                "NEB_Status": r.neb.status.name,
            }
        )

    df = pl.DataFrame(csv_rows)
    df = df.sort(["N_Atoms", "System"])
    df.write_csv(output)

    # 5. Terminal Summary
    console.print("\n[bold]Benchmark Summary[/]")

    # Columns requested: Total Calls, Total Time, Final Convergence, Dimer Travel
    summary_cols = [
        "System",
        "Total_Time_s",
        "Total_Calls",
        "MMF_Status",
        "Dimer_Travel_RMSD",
    ]

    # Filter columns to ensure they exist
    summary_cols = [c for c in summary_cols if c in df.columns]

    pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=console.width)
    console.print(df.select(summary_cols))

    stats_cols = ["Total_Time_s", "Total_Calls", "Dimer_Travel_RMSD"]
    existing = [c for c in stats_cols if c in df.columns]
    summary = df.select([
        *[pl.col(c).mean().alias(f"{c}_mean") for c in existing],
        *[pl.col(c).median().alias(f"{c}_median") for c in existing],
    ])
    pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=console.width)
    console.print(summary)

    console.print(f"\n[bold green]Full detailed metrics saved to: {output}[/]")


if __name__ == "__main__":
    main()

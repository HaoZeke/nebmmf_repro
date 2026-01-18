#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#    "rich",
#    "click",
#    "polars",
#    "ase",
#    "numpy",
# ]
# ///

"""
Optbench NEB Benchmark Tool: Pt Island Diffusion Analysis.

This script compares CINEB and MMF calculations for Platinum island systems.
It performs the following tasks:
1. Parses 'results.dat' for wall-times and saddle point energies.
2. Loads '.con' files to evaluate geometric convergence.
3. Calculates standard RMSD between saddle structures.
4. Exports the compiled dataset to a CSV file for subsequent analysis.
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

# --- Configuration ---

logging.basicConfig(
    level="WARNING",
    format="%(message)s",
    handlers=[RichHandler(show_path=False)],
)

# Regex patterns for EON output files
CALLS_REGEX = re.compile(r"destroyed after (\d+) calls")
EXTREMUM_REGEX = re.compile(r"([\d.-]+)\s+extremum\d+_energy")
TIME_REGEX = re.compile(r"time_seconds\s+([\d.]+)")


class NEBStatus(IntEnum):
    """Maps EON/NEB termination codes to descriptive statuses."""

    GOOD = 0
    INIT = 1
    BAD_MAX_ITERATIONS = 2
    RUNNING = 3
    MAX_UNCERTAINTY = 4
    UNKNOWN = 999

    @classmethod
    def from_code(cls, code: str) -> "NEBStatus":
        try:
            return cls(int(code))
        except (ValueError, TypeError):
            return cls.UNKNOWN


@dataclass
class PathMetrics:
    """Encapsulates energy and geometry metrics for a reaction path."""

    saddle_energy: float | None = None
    saddle_force: float | None = None
    barrier: float | None = None
    saddle_atoms: Atoms | None = None
    path_length: float | None = None
    avg_spacing: float | None = None
    max_spacing: float | None = None


@dataclass
class RunData:
    """Stores the aggregated results for a specific NEB calculation."""

    system: str
    method: str  # 'cineb' or 'mmf'
    time_seconds: float | None = None
    calls: int | None = None
    saddle_energy: float | None = None
    init: PathMetrics | None = None
    final: PathMetrics | None = None
    neb_status: NEBStatus = NEBStatus.UNKNOWN


# --- Core Parsers ---


def get_saddle_energy(content: str) -> float | None:
    """Identifies the maximum energy among all extrema in results.dat."""
    energies = [float(x) for x in EXTREMUM_REGEX.findall(content)]
    return max(energies) if energies else None


def calculate_rmsd(atoms_a: Atoms, atoms_b: Atoms) -> float:
    """Calculates a standard RMSD without rotational or translational alignment."""
    pos_a = atoms_a.get_positions()
    pos_b = atoms_b.get_positions()

    if pos_a.shape != pos_b.shape:
        return -1.0

    diff = pos_a - pos_b
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))


def analyze_path(dat_path: pathlib.Path, con_path: pathlib.Path) -> PathMetrics | None:
    """Extracts energy and geometry statistics from path files."""
    if not dat_path.exists() or not con_path.exists() or con_path.stat().st_size == 0:
        return None

    try:
        # 1. Parse .dat file for energy and force
        lines = dat_path.read_text().splitlines()
        data_map = {}  # index -> (energy, force)

        for line in lines:
            parts = line.strip().split()
            if not parts or not parts[0].isdigit():
                continue
            idx = int(parts[0])
            energy = float(parts[2])
            force = float(parts[3]) if len(parts) > 3 else 0.0
            data_map[idx] = (energy, force)

        if not data_map:
            return None

        # 2. Identify the Saddle Point
        saddle_idx = max(data_map.keys(), key=lambda i: data_map[i][0])
        saddle_energy, saddle_force = data_map[saddle_idx]
        barrier = saddle_energy - data_map[0][0] if 0 in data_map else None

        # 3. Analyze Trajectory Geometry
        traj = ase_read(con_path, index=":")
        if not traj or saddle_idx >= len(traj):
            return None

        saddle_atoms = traj[saddle_idx]
        path_length = 0.0
        spacings = []

        for i in range(len(traj) - 1):
            p1, p2 = traj[i].get_positions(), traj[i + 1].get_positions()
            dist = np.sqrt(np.sum((p2 - p1) ** 2))
            spacings.append(dist)
            path_length += dist

        return PathMetrics(
            saddle_energy=saddle_energy,
            saddle_force=saddle_force,
            barrier=barrier,
            saddle_atoms=saddle_atoms,
            path_length=path_length,
            avg_spacing=np.mean(spacings) if spacings else 0.0,
            max_spacing=np.max(spacings) if spacings else 0.0,
        )
    except Exception:
        return None


def parse_run(path: pathlib.Path) -> RunData | None:
    """Parses a directory containing calculation results."""
    parts = path.parts
    if "cineb" in parts:
        method, sys_idx = "cineb", parts.index("cineb") - 1
    elif "mmf" in parts:
        method, sys_idx = "mmf", parts.index("mmf") - 1
    else:
        return None

    system_name = parts[sys_idx]
    results_file = path / "results.dat"

    # Analyze path files
    init_metrics = analyze_path(path / "neb_000.dat", path / "neb_path_000.con")
    final_metrics = analyze_path(path / "neb.dat", path / "neb.con")

    time_val, s_energy, status = None, None, NEBStatus.UNKNOWN
    if results_file.exists():
        txt = results_file.read_text(errors="ignore")
        t_match = TIME_REGEX.search(txt)
        time_val = float(t_match.group(1)) if t_match else None
        s_energy = get_saddle_energy(txt)

        first_line = txt.splitlines()[0] if txt else ""
        if "termination_reason" in first_line:
            status = NEBStatus.from_code(first_line.split()[0])

    calls = None
    for log in path.glob("*.log"):
        l_txt = log.read_text(errors="ignore")
        c_match = CALLS_REGEX.search(l_txt)
        if c_match:
            calls = int(c_match.group(1))
            break

    return RunData(
        system_name,
        method,
        time_val,
        calls,
        s_energy,
        init_metrics,
        final_metrics,
        status,
    )


# --- Main Logic ---


@click.command()
@click.argument(
    "search_path", type=click.Path(exists=True, path_type=pathlib.Path), default="."
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=pathlib.Path),
    default="optbench_results.csv",
)
def main(search_path: pathlib.Path, output: pathlib.Path):
    """Crawls directories for Pt island NEB runs and generates a comparison CSV."""
    console = Console()
    results_files = list(search_path.rglob("results.dat"))
    run_dirs = [p.parent for p in results_files]

    if not run_dirs:
        console.print("[red]No calculation results found.[/]")
        return

    console.print(f"[bold green]Parsing {len(run_dirs)} runs...[/]")
    runs = []
    with ThreadPoolExecutor() as exc:
        futures = [exc.submit(parse_run, d) for d in run_dirs]
        for f in futures:
            res = f.result()
            if res:
                runs.append(res)

    data_dicts = []
    for r in runs:
        data_dicts.append(
            {
                "system": r.system,
                "method": r.method,
                "time": r.time_seconds,
                "calls": r.calls,
                "energy": r.saddle_energy,
                "termination": r.neb_status.name,
                "init": r.init,
                "final": r.final,
            }
        )

    df = pl.DataFrame(
        data_dicts, schema_overrides={"init": pl.Object, "final": pl.Object}
    )
    cineb = df.filter(pl.col("method") == "cineb").drop("method")
    mmf = df.filter(pl.col("method") == "mmf").drop("method")
    joined = cineb.join(mmf, on="system", suffix="_mmf", how="outer")

    final_rows = []
    with Progress(console=console) as progress:
        task = progress.add_task("[cyan]Processing Data...", total=len(joined))

        for row in joined.iter_rows(named=True):
            sys_name = row["system"]

            rmsd_saddle = None
            if row["final"] and row["final_mmf"]:
                rmsd_saddle = calculate_rmsd(
                    row["final"].saddle_atoms, row["final_mmf"].saddle_atoms
                )

            final_rows.append(
                {
                    "System": sys_name,
                    "E_Diff": (row["energy"] - row["energy_mmf"])
                    if row["energy"] and row["energy_mmf"]
                    else None,
                    "RMSD_Saddle": rmsd_saddle,
                    "Ratio_Calls": (row["calls"] / row["calls_mmf"])
                    if row["calls"] and row["calls_mmf"]
                    else None,
                    "Time_CINEB": row["time"],
                    "Time_MMF": row["time_mmf"],
                    "Calls_CINEB": row["calls"],
                    "Calls_MMF": row["calls_mmf"],
                    "Force_CINEB": row["final"].saddle_force if row["final"] else None,
                    "Force_MMF": row["final_mmf"].saddle_force
                    if row["final_mmf"]
                    else None,
                }
            )
            progress.advance(task)

    result_df = pl.DataFrame(final_rows).sort("System")
    result_df.write_csv(output)

    pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=console.width)
    console.print(result_df)
    console.print(f"[bold green]Analysis saved to {output}[/]")


if __name__ == "__main__":
    main()

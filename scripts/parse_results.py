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
Advanced NEB Benchmark Tool: Geometry & Energy Analysis.

Compare CINEB vs MMF runs by:
1. Parsing 'results.dat' for wall-times and saddle energies.
2. Loading 'sp.con' files to compare geometries.
3. Using 'rgpycrumbs' to robustly align structures (IRA/Procrustes) and calculate RMSD.
4. Exporting a rich dataset for R analysis.
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
    # Fallback if running outside the specific env, primarily for the linter/demo
    logging.warning(
        "Could not import 'rgpycrumbs'. Geometric alignment will be skipped."
    )
    align_structure_robust = None

# --- Configuration ---

logging.basicConfig(
    level="WARNING",  # Keep logs quiet, focus on the table/progress
    format="%(message)s",
    handlers=[RichHandler(show_path=False)],
)

# Regex patterns
CALLS_REGEX = re.compile(r"destroyed after (\d+) calls")
EXTREMUM_REGEX = re.compile(r"([\d.-]+)\s+extremum\d+_energy")
TIME_REGEX = re.compile(r"time_seconds\s+([\d.]+)")


class NEBStatus(IntEnum):
    """Maps EON/NEB termination codes to human-readable statuses."""

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
    """Encapsulates physics (from .dat) and geometry (from .con) metrics."""

    # From .dat
    saddle_energy: float | None = None
    saddle_force: float | None = None
    barrier: float | None = None

    # From .con
    saddle_atoms: Atoms | None = None
    path_length: float | None = None
    avg_spacing: float | None = None
    max_spacing: float | None = None


@dataclass
class RunData:
    system: str
    method: str  # 'cineb' or 'mmf'
    time_seconds: float | None = None
    calls: int | None = None
    saddle_energy: float | None = None
    init: PathMetrics | None = None
    final: PathMetrics | None = None
    neb_status: NEBStatus = NEBStatus.UNKNOWN


# --- Reaction Metadata Mapping ---

# REACTION_MAP = {
#     "01": "HCN → HNC",
#     "02": "HCCH → CCH₂",
#     "03": "H₂CO → H₂ + CO",
#     "04": "CH₃O → CH₂OH",
#     "05": "cyclopropyl ring opening",
#     "06": "bicyclo[1.1.0]butane → trans-butadiene",
#     "08": "formyloxyethyl 1,2-migration",
#     "09": "parent Diels-Alder cycloaddition",
#     "10": "s-tetrazine → 2HCN + N₂",
#     "11": "trans-butadiene → cis-butadiene",
#     "12": "CH₃CH₃ → CH₂CH₂ + H₂",
#     "13": "CH₃CH₂F → CH₂CH₂ + HF",
#     "14": "acetaldehyde keto-enol tautomerism",
#     "15": "HOCl → HCl + CO",
#     "16": "H₂O + PO₃⁻ → H₂PO₄⁻",
#     "17": "CH₂CHCH₂CH₂CHO Claisen rearrangement",
#     "18": "SiH₂ + CH₃CH₃ → SiH₃CH₂CH₃",
#     "19": "HNCCS → HNC + CS",
#     "20": "HCONH₃⁺ → NH₄⁺ + CO",
#     "21": "acrolein rotational TS",
#     "22": "HCONHOH → HCOHNHO",
#     "23": "HNC + H₂ → H₂CNH",
#     "24": "H₂CNH → HCNH₂",
#     "25": "HCNH₂ → HCN + H₂",
# }

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

# --- Core Parsers ---


def get_saddle_energy(content: str) -> float | None:
    """Finds the highest energy among all extrema in results.dat (the Saddle Point)."""
    energies = [float(x) for x in EXTREMUM_REGEX.findall(content)]
    return max(energies) if energies else None


def parse_run(path: pathlib.Path) -> RunData | None:
    """Parses a single run directory for physics and geometry data."""
    parts = path.parts
    # Heuristic to find system/method from path: .../system_name/method_name
    if "cineb" in parts:
        method = "cineb"
        sys_idx = parts.index("cineb") - 1
    elif "mmf" in parts:
        method = "mmf"
        sys_idx = parts.index("mmf") - 1
    else:
        return None

    system_name = parts[sys_idx]

    # 1. Parse results.dat
    time_val, s_energy = None, None
    term_reason = "Unknown"
    results_file = path / "results.dat"
    # 2. Get Initial Path Metrics (neb_000)
    init_metrics = analyze_path(path / "neb_000.dat", path / "neb_path_000.con")

    # 3. Get Final Path Metrics (neb)
    final_metrics = analyze_path(path / "neb.dat", path / "neb.con")

    if results_file.exists():
        txt = results_file.read_text(errors="ignore")

        # Get Time
        t_match = TIME_REGEX.search(txt)
        if t_match:
            time_val = float(t_match.group(1))

        # Get Saddle Energy
        s_energy = get_saddle_energy(txt)

        # Get Termination (first line typically)
        first_line = txt.splitlines()[0] if txt else ""
        if "termination_reason" in first_line:
            term_reason = NEBStatus.from_code(first_line.split()[0])

    # 2. Parse Calls from logs
    calls = None
    for log in path.glob("*.log"):
        try:
            l_txt = log.read_text(errors="ignore")
            c_match = CALLS_REGEX.search(l_txt)
            if c_match:
                calls = int(c_match.group(1))
                break
        except:
            pass

    # 3. Load Geometry (Saddle Point)
    atoms = None
    sp_file = path / "sp.con"
    if sp_file.exists():
        try:
            # ASE usually handles .con if configured, otherwise assumed generic format
            # If standard ASE fails, we might need 'format="eon"' if installed
            atoms = ase_read(sp_file)
            assert atoms == final_metrics.saddle
        except Exception:
            pass

    return RunData(
        system_name,
        method,
        time_val,
        calls,
        s_energy,
        init_metrics,
        final_metrics,
        term_reason,
    )


def calculate_rmsd(atoms_a: Atoms, atoms_b: Atoms) -> float:
    """Calculates RMSD using the user's robust alignment (IRA/Procrustes)."""
    if not align_structure_robust:
        return -1.0

    # We clone to avoid modifying the cached objects in place if we need them later
    mobile = atoms_b.copy()
    ref = atoms_a.copy()

    config = IRAConfig(enabled=True, kmax=1.8)

    # This aligns 'mobile' onto 'ref' in-place
    align_structure_robust(ref, mobile, config)

    # Calculate RMSD
    pos_a = ref.get_positions()
    pos_b = mobile.get_positions()
    diff = pos_a - pos_b
    rmsd = np.sqrt((diff * diff).sum() / len(ref))
    return rmsd


def analyze_path(dat_path: pathlib.Path, con_path: pathlib.Path) -> PathMetrics | None:
    """
    Parses .dat for Energy/Force and .con for Geometry stats.
    """
    if not dat_path.exists() or not con_path.exists() or con_path.stat().st_size == 0:
        return None

    try:
        # 1. Parse .dat file (img, rxn_coord, energy, force)
        lines = dat_path.read_text().splitlines()
        data_map = {}  # index -> (energy, force)

        for line in lines:
            parts = line.strip().split()
            if not parts or not parts[0].isdigit():
                continue
            try:
                idx = int(parts[0])
                energy = float(parts[2])
                # Column 3 is usually f_para or force magnitude
                # Handle cases where column 3 might be missing (though rare in EON)
                force = float(parts[3]) if len(parts) > 3 else 0.0
                data_map[idx] = (energy, force)
            except (ValueError, IndexError):
                continue

        if not data_map:
            return None

        # 2. Identify Saddle (Max Energy Image)
        saddle_idx = max(data_map.keys(), key=lambda i: data_map[i][0])
        saddle_energy, saddle_force = data_map[saddle_idx]

        # Barrier (Saddle - Reactant)
        barrier = None
        if 0 in data_map:
            barrier = saddle_energy - data_map[0][0]

        # 3. Load Trajectory for Spacing/Length
        # We need the atoms to calculate geometric path length and find the saddle structure
        traj = ase_read(con_path, index=":")
        if not traj or saddle_idx >= len(traj):
            return None

        saddle_atoms = traj[saddle_idx]

        # 4. Calculate Path Statistics (Geometry)
        path_length = 0.0
        spacings = []

        for i in range(len(traj) - 1):
            p1 = traj[i].get_positions()
            p2 = traj[i + 1].get_positions()
            # Standard NEB distance (Euclidean norm of flattened coordinate vector)
            d = np.sqrt(np.sum((p2 - p1) ** 2))
            spacings.append(d)
            path_length += d

        avg_spacing = np.mean(spacings) if spacings else 0.0
        max_spacing = np.max(spacings) if spacings else 0.0

        return PathMetrics(
            saddle_energy=saddle_energy,
            saddle_force=saddle_force,
            barrier=barrier,
            saddle_atoms=saddle_atoms,
            path_length=path_length,
            avg_spacing=avg_spacing,
            max_spacing=max_spacing,
        )

    except Exception:
        return None


# --- Main Logic ---


@click.command()
@click.argument(
    "search_path", type=click.Path(exists=True, path_type=pathlib.Path), default="."
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=pathlib.Path),
    default="benchmark_results.csv",
    help="Output CSV path.",
)
def main(search_path: pathlib.Path, output: pathlib.Path):
    """
    Crawls SEARCH_PATH for matched CINEB/MMF runs, aligns geometries, and outputs a CSV.
    """
    console = Console()

    # 1. Find Directories
    results_files = list(search_path.rglob("results.dat"))
    run_dirs = [p.parent for p in results_files]

    if not run_dirs:
        console.print("[red]No results.dat files found![/]")
        return

    # 2. Parse Data in Parallel
    console.print(f"[bold green]Parsing {len(run_dirs)} runs...[/]")
    runs = []
    with ThreadPoolExecutor() as exc:
        # We process files to get raw data
        futures = [exc.submit(parse_run, d) for d in run_dirs]
        for f in futures:
            res = f.result()
            if res:
                runs.append(res)

    # 3. Aggregate into Polars
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

    # 4. Pivot and Compare
    # Split into CINEB and MMF frames
    cineb = df.filter(pl.col("method") == "cineb").drop("method")
    mmf = df.filter(pl.col("method") == "mmf").drop("method")

    # Join on System
    joined = cineb.join(mmf, on="system", suffix="_mmf", how="outer")

    # 5. Calculate Physics & Geometry Deltas
    # We iterate over the joined rows to do the heavy geometric lifting (RMSD)

    final_rows = []

    with Progress(console=console) as progress:
        task = progress.add_task(
            "[cyan]Aligning Structures & Calculating Stats...", total=len(joined)
        )

        for row in joined.iter_rows(named=True):
            sys_name = row["system"]
            label = REACTION_MAP.get(sys_name.split("_")[0], "Unknown Reaction")

            # --- Chemical Metadata (from CINEB atoms, or MMF if CINEB missing) ---
            ref_atoms = row["init"].saddle_atoms or row["init_mmf"].saddle_atoms
            if ref_atoms:
                natoms = len(ref_atoms)
                formula = ref_atoms.get_chemical_formula()
                mass = sum(ref_atoms.get_masses())
            else:
                natoms, formula, mass = 0, "Unknown", 0.0

            # --- Calculate "Distance Traveled" (Initial Guess -> Final Saddle) ---
            # This quantifies how "bad" the initial guess was.
            rmsd_travel_cineb = None
            if row["init"].saddle_atoms and row["final"].saddle_atoms:
                try:
                    rmsd_travel_cineb = calculate_rmsd(
                        row["init"].saddle_atoms, row["final"].saddle_atoms
                    )
                except Exception:
                    pass
            rmsd_travel_mmf = None
            if row["init_mmf"].saddle_atoms and row["final_mmf"].saddle_atoms:
                try:
                    rmsd_travel_mmf = calculate_rmsd(
                        row["init_mmf"].saddle_atoms, row["final_mmf"].saddle_atoms
                    )
                except Exception:
                    pass

            f_init_cineb = row["init"].saddle_force if row["init"] else None
            f_init_mmf = row["init_mmf"].saddle_force if row["init_mmf"] else None
            f_final_cineb = row["final"].saddle_force if row["final"] else None
            f_final_mmf = row["final_mmf"].saddle_force if row["final_mmf"] else None

            # --- Geometric Comparison ---
            rmsd_finals = None
            if row["final"].saddle_atoms and row["final_mmf"].saddle_atoms:
                try:
                    # Align MMF (mobile) to CINEB (ref)
                    rmsd_finals = calculate_rmsd(
                        row["final"].saddle_atoms, row["final_mmf"].saddle_atoms
                    )
                except Exception as e:
                    logging.warning(f"Alignment failed for {sys_name}: {e}")

            # --- Energy Difference ---
            e_diff = None
            if row["energy"] is not None and row["energy_mmf"] is not None:
                e_diff = row["energy"] - row["energy_mmf"]

            # --- Barrier Difference ---
            # Using parsed PathMetrics barrier if available
            b_cineb = row["final"].barrier
            b_mmf = row["final_mmf"].barrier
            b_diff = None
            if b_cineb is not None and b_mmf is not None:
                b_diff = b_cineb - b_mmf

            # --- Ratios & Efficiency ---
            ratio_calls = None
            if row["calls"] and row["calls_mmf"]:
                ratio_calls = row["calls"] / row["calls_mmf"]

            ratio_time = None
            if row["time"] and row["time_mmf"]:
                ratio_time = row["time"] / row["time_mmf"]

            # --- Construct Final Row ---
            final_rows.append(
                {
                    "System": sys_name,
                    "Reaction": label,
                    "Formula": formula,
                    "N_Atoms": natoms,
                    "Total_Mass": mass,
                    # Comparison
                    "E_Diff": e_diff,
                    "Barrier_Diff": b_diff,
                    "RMSD_Saddle": rmsd_finals,
                    "Ratio_Calls": ratio_calls,
                    "Ratio_Time": ratio_time,
                    "RMSD_Init_Final_CINEB": rmsd_travel_cineb,
                    "RMSD_Init_Final_MMF": rmsd_travel_mmf,
                    # CINEB Specifics
                    "Barrier_CINEB": b_cineb,
                    "Force_CINEB": row["final"].saddle_force,
                    "Path_Len_CINEB": row["final"].path_length,
                    "Calls_CINEB": row["calls"],
                    "Time_CINEB": row["time"],
                    "Term_CINEB": row["termination"],
                    "Force_Init_CINEB": f_init_cineb,
                    "Force_Final_CINEB": f_final_cineb,
                    # MMF Specifics
                    "Barrier_MMF": b_mmf,
                    "Force_MMF": row["final_mmf"].saddle_force,
                    "Calls_MMF": row["calls_mmf"],
                    "Time_MMF": row["time_mmf"],
                    "Term_MMF": row["termination_mmf"],
                    "Force_Init_MMF": f_init_mmf,
                    "Force_Final_MMF": f_final_mmf,
                }
            )
            progress.advance(task)

    # 6. Create Final DataFrame & Save
    result_df = pl.DataFrame(final_rows)

    # Sort by complexity (N_Atoms) then System name
    result_df = result_df.sort(["N_Atoms", "System"])

    # Write to CSV
    result_df.write_csv(output)

    # Print a quick summary table
    pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=console.width)
    console.print(
        result_df.sort(by="System").select(
            [
                "System",
                "Formula",
                "Time_CINEB",
                "Time_MMF",
                "E_Diff",
                "RMSD_Saddle",
                "Calls_CINEB",
                "Calls_MMF",
            ]
        )
    )
    console.print(f"[bold green]Saved detailed analysis to {output}[/]")


if __name__ == "__main__":
    main()

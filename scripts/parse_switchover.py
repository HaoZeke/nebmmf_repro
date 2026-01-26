#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "rich",
#   "click",
#   "polars",
# ]
# ///

"""
NEB-MMF Switch Analyzer.

Parses EON client logs to collate NEB<->MMF switching statistics.
Extracts specific physical reasons (Force/Alignment thresholds) for
fallback events to allow statistical analysis of MMF stability.
"""

import re
import pathlib
import click
import polars as pl
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, SpinnerColumn, TextColumn
import logging

# --- Configuration ---

logging.basicConfig(
    level="INFO",
    format="%(message)s",
    handlers=[RichHandler(show_path=False)],
)

log = logging.getLogger("SwitchParser")

# --- Regex Patterns ---

# Start of a run block
RX_START = re.compile(r"^EON Client")
RX_BUILD = re.compile(r"BUILD DATE:\s+(.*)")
RX_DIR = re.compile(r"DIR:\s+(.*)")

# Switch Logic: NEB -> MMF
# "Triggering MMF.  Force: 1.2333, Threshold: 2.4131 (0.50x baseline)"
RX_TRIGGER = re.compile(
    r"Triggering MMF\.\s+Force:\s+(?P<force>[\d\.]+),\s+Threshold:\s+(?P<threshold>[\d\.]+)"
)

# Switch Logic: MMF -> NEB (The Backoff)
# "MMF backoff (status=0). Force: 1.2333 -> 1.4775, Alignment:  0.923."
RX_BACKOFF = re.compile(
    r"MMF backoff \(status=(?P<status>\d+)\)\.\s+Force:\s+(?P<f_old>[\d\.]+)\s+->\s+(?P<f_new>[\d\.]+),\s+Alignment:\s+(?P<align>[\d\.]+)(?=\.)"
)

# Convergence
RX_CONVERGED_MMF = re.compile(r"NEB converged after MMF")
RX_CONVERGED_NEB = re.compile(r"NEB converged")


# --- Parsing Logic ---


def parse_log_stream(content: str) -> pl.DataFrame:
    """
    Parses concatenated log content into a structured Polars DataFrame.
    One row per Run, with nested lists for switch events.
    """
    runs = []

    # State tracking
    current_run = None
    lines = content.splitlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue

        # Detect Start of New Run
        if RX_START.search(line):
            if current_run:
                runs.append(current_run)

            # Reset
            current_run = {
                "system": "Unknown",
                "build_date": None,
                "n_triggers": 0,
                "n_backoffs": 0,
                "trigger_forces": [],  # List[Float]
                "backoff_forces_pre": [],  # List[Float]
                "backoff_forces_post": [],  # List[Float]
                "backoff_alignments": [],  # List[Float]
                "final_state": "Incomplete",
            }
            continue

        # If we haven't found a 'start' tag yet, skip metadata parsing
        if current_run is None:
            continue

        # Metadata
        if m := RX_BUILD.search(line):
            current_run["build_date"] = m.group(1)

        elif m := RX_DIR.search(line):
            # Extract relevant path segment (e.g., '06_bicyclobutane')
            path_str = m.group(1)
            parts = path_str.split("/")
            # Heuristic: Grab the segment before 'mmf' or 'cineb' if possible
            if "mmf" in parts:
                idx = parts.index("mmf")
                current_run["system"] = parts[idx - 1] if idx > 0 else path_str
            else:
                current_run["system"] = parts[-1]

        # Trigger Event (NEB -> MMF)
        elif m := RX_TRIGGER.search(line):
            current_run["n_triggers"] += 1
            current_run["trigger_forces"].append(float(m.group("force")))

        # Backoff Event (MMF -> NEB)
        elif m := RX_BACKOFF.search(line):
            current_run["n_backoffs"] += 1
            current_run["backoff_forces_pre"].append(float(m.group("f_old")))
            current_run["backoff_forces_post"].append(float(m.group("f_new")))
            current_run["backoff_alignments"].append(float(m.group("align")))

        # Final Status
        elif RX_CONVERGED_MMF.search(line):
            current_run["final_state"] = "Converged (MMF)"
        elif RX_CONVERGED_NEB.search(line):
            current_run["final_state"] = "Converged (NEB)"

    # Append the last run
    if current_run:
        runs.append(current_run)

    # Convert to Polars
    # We must explicitely tell Polars that the lists contain Floats,
    # otherwise empty lists might be inferred as List[Null] or similar.
    schema = {
        "system": pl.String,
        "build_date": pl.String,
        "n_triggers": pl.UInt32,
        "n_backoffs": pl.UInt32,
        "trigger_forces": pl.List(pl.Float64),
        "backoff_forces_pre": pl.List(pl.Float64),
        "backoff_forces_post": pl.List(pl.Float64),
        "backoff_alignments": pl.List(pl.Float64),
        "final_state": pl.String,
    }

    return pl.DataFrame(runs, schema=schema)


# --- CLI ---


@click.command()
@click.argument("root_path", type=click.Path(exists=True, path_type=pathlib.Path))
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=pathlib.Path),
    default="switches.parquet",
    help="Output file path (parquet/csv/json).",
)
@click.option(
    "--csv", is_flag=True, help="Force CSV output (flattens lists to strings)."
)
def main(root_path: pathlib.Path, output: pathlib.Path, csv: bool):
    """
    Analyzes EON logs for NEB/MMF switching behavior.

    ROOT_PATH can be a single log file or a directory containing *.log files.
    """
    console = Console()

    # Gather
    console.print(f"[bold]Scanning {root_path} for MMF logs...[/]")

    found_files = []
    # Find all 'mmf' directories first
    mmf_dirs = [p for p in root_path.rglob("mmf") if p.is_dir()]

    if not mmf_dirs:
        # Fallback: maybe the user pointed directly at a log file?
        if root_path.is_file():
            found_files = [root_path]
        else:
            console.print("[red]No 'mmf' subdirectories found![/]")
            return
    else:
        for d in mmf_dirs:
            # Grab all .log files in this mmf folder
            found_files.extend(list(d.glob("*.log")))

    if not found_files:
        console.print("[red]Found 'mmf' folders but no .log files inside them.[/]")
        return

    full_content = []

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        task = progress.add_task(
            f"Reading {len(found_files)} logs...", total=len(found_files)
        )

        for f in found_files:
            try:
                full_content.append(f.read_text(errors="ignore"))
            except Exception as e:
                log.warning(f"Could not read {f}: {e}")
            progress.advance(task)

    # Parse
    console.print("[bold blue]Parsing event streams...[/]")
    combined_text = "\n".join(full_content)
    df = parse_log_stream(combined_text)

    # Output
    if df.is_empty():
        console.print("[yellow]No runs detected in logs.[/]")
        return

    # Sort by complexity (triggers)
    df = df.sort("n_backoffs", descending=True)

    # Pretty print summary to console
    console.print("\n[bold]Run Summary[/]")
    with pl.Config(tbl_rows=-1, tbl_cols=-1):
        df = df.sort("system")
        console.print(df.select(["system", "n_triggers", "n_backoffs", "final_state"]))

    # Save
    if csv or output.suffix == ".csv":
        # Polars writes lists as strings in CSV "[1.0, 2.0]", which is usually fine for text analysis
        # but Parquet is better for numerical analysis in Python/R.
        df.write_csv(output)
        console.print(f"\n[green]Saved CSV to {output}[/]")
    else:
        df.write_parquet(output)
        console.print(f"\n[green]Saved Parquet to {output}[/]")

    # Hint for user
    if df["n_backoffs"].sum() > 0:
        console.print(
            "\n[italic dim]Tip: Use 'df.explode(\"backoff_alignments\")' to plot a histogram of failure angles.[/]"
        )


if __name__ == "__main__":
    main()

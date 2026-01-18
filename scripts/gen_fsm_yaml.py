#!/usr/bin/env python3

# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "click",
#     "rich",
#     "pyyaml"
# ]
# ///
#
# Usage:
# uv run gen_fsm_yaml.py -i resources/icFSM/baker -o baker_config.yaml

from collections import defaultdict
from pathlib import Path

import click
import yaml
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

console = Console()


def discover_systems(base_dir: Path) -> dict:
    """
    Recursively scans the directory for system patterns like 01_hcn.
    Looks for initial/ipath_*.con and ts.xyz files.
    """
    systems_data = defaultdict(dict)

    # Pattern to match directories like '01_hcn' or 'ethanal_03'
    # We look for folders that contain subfolders or specific NEB files
    system_dirs = [d for d in base_dir.iterdir() if d.is_dir()]

    for sdir in system_dirs:
        system_id = sdir.name

        # 1. Search for Reactant/Product in 'initial' subfolder
        initial_path = sdir / "initial"
        if initial_path.exists():
            con_files = sorted(list(initial_path.glob("ipath_*.con")))
            if len(con_files) >= 2:
                # Assuming standard 000 is reactant and last is product
                # or 001/002 based on your previous config
                systems_data[system_id]["reactant"] = str(con_files[0].resolve())
                systems_data[system_id]["product"] = str(con_files[-1].resolve())
            elif len(con_files) == 1:
                systems_data[system_id]["reactant"] = str(con_files[0].resolve())

        # 2. Search for Saddle (ts.xyz)
        saddle_file = sdir / "ts.xyz"
        if saddle_file.exists():
            systems_data[system_id]["saddle"] = str(saddle_file.resolve())

        # 3. Add default flags
        # if system_id in systems_data:
        # systems_data[system_id]["use_ira"] = False

    return systems_data


def write_config(data: dict, output_file: str):
    """Outputs the gathered data to a YAML file."""
    final_structure = {"systems": dict(data)}
    try:
        with open(output_file, "w") as f:
            yaml.dump(final_structure, f, sort_keys=False, default_flow_style=False)

        console.print(
            Panel(
                f"[green]Generated config for {len(data)} systems.[/green]\nFile: [yellow]{output_file}[/yellow]"
            )
        )

        # Summary Table
        table = Table(title="Discovery Summary")
        table.add_column("System", style="cyan")
        table.add_column("Reactant")
        table.add_column("Product")
        table.add_column("Saddle")

        for name, paths in data.items():
            table.add_row(
                name,
                "[green]ok[/green]" if "reactant" in paths else "[red]-[/red]",
                "[green]ok[/green]" if "product" in paths else "[red]-[/red]",
                "[green]ok[/green]" if "saddle" in paths else "[red]-[/red]",
            )
        console.print(table)

    except Exception as e:
        console.print(f"[bold red]Error writing YAML:[/bold red] {e}")


@click.command()
@click.option(
    "--input-dir", "-i", type=click.Path(exists=True, path_type=Path), required=True
)
@click.option("--output-file", "-o", default="neb_systems.yaml")
def main(input_dir: Path, output_file: str):
    """Scans icFSM/mlFSM style directories to generate Snakemake system configs."""
    console.rule("[bold magenta]NEB System Discovery[/bold magenta]")

    systems = discover_systems(input_dir)
    if not systems:
        console.print(
            "[yellow]No systems found matching the directory pattern.[/yellow]"
        )
        return

    write_config(systems, output_file)


if __name__ == "__main__":
    main()

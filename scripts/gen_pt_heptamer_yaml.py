#!/usr/bin/env python3

# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "click",
#     "rich",
#     "pyyaml"
# ]
# ///
# Usage:
# uv run gen_pt_heptamer_yaml.py -i resources/pt-island-con -o pthept_config.yaml

import re
from pathlib import Path

import click
import yaml
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

console = Console()


def get_product_index(path: Path) -> int:
    """Extracts the numerical index from a product_N.con filename."""
    match = re.search(r"product_(\d+)\.con", path.name)
    return int(match.group(1)) if match else -1


def discover_systems(base_dir: Path) -> dict:
    """
    Scans the directory for a reactant.con and multiple product_*.con files.
    Maps them into a systems dictionary compatible with the NEB workflow.
    """
    systems_data = {}

    # 1. Locate the common reactant
    reactant_file = base_dir / "reactant.con"
    if not reactant_file.exists():
        console.print(f"[bold red]Error:[/bold red] Could not find {reactant_file}")
        return {}

    # 2. Find all product files
    product_files = sorted(list(base_dir.glob("product_*.con")))
    product_files.sort(key=get_product_index)

    for p_file in product_files:
        # Extract the index from 'product_N.con' to create a clean system ID
        match = re.search(r"product_(\d+)\.con", p_file.name)
        if match:
            idx = int(match.group(1))
            # Format system ID as pt_prod_00, pt_prod_01, etc.
            system_id = f"pt_prod_{idx:02d}"

            systems_data[system_id] = {
                "reactant": str(reactant_file.relative_to(base_dir.parent.parent)),
                "product": str(p_file.relative_to(base_dir.parent.parent)),
            }

    return systems_data


def write_config(data: dict, output_file: str):
    """Outputs the gathered data to a YAML file."""
    final_structure = {"systems": dict(data)}
    try:
        with open(output_file, "w") as f:
            yaml.dump(final_structure, f, sort_keys=False, default_flow_style=False)

        console.print(
            Panel(
                f"[green]Generated config for {len(data)} systems.[/green]\n"
                f"File: [yellow]{output_file}[/yellow]"
            )
        )

        # Summary Table
        table = Table(title="Discovery Summary")
        table.add_column("System ID", style="cyan")
        table.add_column("Product File", style="green")
        table.add_column("IRA Enabled", style="magenta")

        for name, entry in data.items():
            product_filename = Path(entry["product"]).name
            table.add_row(name, product_filename)
        console.print(table)

    except Exception as e:
        console.print(f"[bold red]Error writing YAML:[/bold red] {e}")


@click.command()
@click.option(
    "--input-dir",
    "-i",
    type=click.Path(exists=True, path_type=Path),
    default="resources/pt-island-con",
    help="Directory containing reactant.con and product_*.con files.",
)
@click.option("--output-file", "-o", default="config/systems_config.yaml")
def main(input_dir: Path, output_file: str):
    """Generates Snakemake system configs from Pt-island style directories."""
    console.rule("[bold magenta]Pt-Island System Discovery[/bold magenta]")

    systems = discover_systems(input_dir)
    if not systems:
        console.print(
            "[yellow]No product files found in the specified directory.[/yellow]"
        )
        return

    write_config(systems, output_file)


if __name__ == "__main__":
    main()

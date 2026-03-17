#!/usr/bin/env python3
"""Generate all supplementary material figures for the OCI-NEB paper revision.

Produces 2D RMSD landscape plots for each Baker system (both CI-NEB and MMF
methods), and optionally runs the speedup violin plot script.
"""
import argparse
import subprocess
import sys
from pathlib import Path

BAKER_SYSTEMS = [
    "01_hcn",
    "02_hcch",
    "03_h2co",
    "04_ch3o",
    "05_cyclopropyl",
    "06_bicyclobutane",
    # 07 is skipped
    "08_formyloxyethyl",
    "09_parentdielsalder",
    "10_tetrazine",
    "11_trans_butadiene",
    "12_ethane_h2_abstraction",
    "13_hf_abstraction",
    "14_vinyl_alcohol",
    "15_hocl",
    "16_h2po4_anion",
    "17_claisen",
    "18_silylene_insertion",
    "19_hnccs",
    "20_hconh3_cation",
    "21_acrolein_rot",
    "22_hconhoh",
    "23_hcn_h2",
    "24_h2cnh",
    "25_hcnh2",
]

METHODS = ["cineb", "mmf"]

DEFAULT_RESULTS_DIR = "eonRuns/results/03_neb"
DEFAULT_OUTPUT_DIR = (
    "/home/rgoswami/Git/Github/TeX/roneb_tex/v2"
    "/FrontiersChemistry/imgs/suppl/baker_2d"
)


def build_landscape_cmd(
    con_file: Path,
    output_file: Path,
    system_dir: Path,
) -> list[str]:
    """Build the rgpycrumbs CLI command for a single landscape plot."""
    return [
        sys.executable,
        "-m",
        "rgpycrumbs.cli",
        "eon",
        "plt-neb",
        "--con-file",
        str(con_file),
        "--output-file",
        str(output_file),
        "--plot-type",
        "landscape",
        "--rc-mode",
        "rmsd",
        "--surface-type",
        "grad_imq",
        "--show-pts",
        "--landscape-path",
        "all",
        "--input-dat-pattern",
        str(system_dir / "neb_*.dat"),
        "--input-path-pattern",
        str(system_dir / "neb_path*.con"),
        "--ira-kmax",
        "14",
        "--plot-structures",
        "crit_points",
        "--dpi",
        "300",
        "--show-legend",
    ]


def generate_landscape_plots(
    results_dir: Path,
    output_dir: Path,
) -> list[Path]:
    """Generate 2D RMSD landscape plots for all Baker systems and methods."""
    generated = []
    failures = []

    for system in BAKER_SYSTEMS:
        for method in METHODS:
            system_dir = results_dir / system / method
            con_file = system_dir / "neb.con"

            if not system_dir.is_dir():
                print(f"SKIP: {system_dir} does not exist")
                continue
            if not con_file.is_file():
                print(f"SKIP: {con_file} not found")
                continue

            out_subdir = output_dir / system
            out_subdir.mkdir(parents=True, exist_ok=True)
            output_file = out_subdir / f"{method}_landscape.png"

            cmd = build_landscape_cmd(con_file, output_file, system_dir)
            print(f"  {system}/{method} -> {output_file}")

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
            )

            if result.returncode != 0:
                failures.append((system, method, result.stderr.strip()))
                print(f"    FAILED (exit {result.returncode})")
                if result.stderr:
                    # Show last few lines of stderr for diagnostics
                    for line in result.stderr.strip().splitlines()[-5:]:
                        print(f"      {line}")
            else:
                generated.append(output_file)
                print("    OK")

    if failures:
        print(f"\n{len(failures)} plot(s) failed:")
        for system, method, err in failures:
            short = err.splitlines()[-1] if err else "(no stderr)"
            print(f"  {system}/{method}: {short}")

    return generated


def run_speedup_violin(output_dir: Path) -> Path | None:
    """Run the speedup violin plot script if it exists."""
    # The script lives in the v2 suppl images directory
    script = output_dir.parent / "make_speedup_violin.py"
    if not script.is_file():
        print(f"SKIP: speedup violin script not found at {script}")
        return None

    print(f"Running speedup violin script: {script}")
    result = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(f"  FAILED (exit {result.returncode})")
        if result.stderr:
            for line in result.stderr.strip().splitlines()[-5:]:
                print(f"    {line}")
        return None

    if result.stdout:
        print(result.stdout.strip())

    violin_out = output_dir.parent / "speedup_ratio_violin.png"
    if violin_out.is_file():
        return violin_out
    return None


def print_manifest(landscape_files: list[Path], violin_file: Path | None) -> None:
    """Print a manifest of all generated files."""
    all_files = list(landscape_files)
    if violin_file is not None:
        all_files.append(violin_file)

    print("\n" + "=" * 60)
    print(f"MANIFEST: {len(all_files)} file(s) generated")
    print("=" * 60)
    for f in sorted(all_files):
        print(f"  {f}")
    print("=" * 60)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate supplementary figures for the OCI-NEB paper.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path(DEFAULT_RESULTS_DIR),
        help="Path to NEB results (default: %(default)s)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(DEFAULT_OUTPUT_DIR),
        help="Output directory for 2D landscape plots (default: %(default)s)",
    )
    parser.add_argument(
        "--skip-violin",
        action="store_true",
        help="Skip the speedup violin plot generation",
    )
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    output_dir = args.output_dir.resolve()

    if not results_dir.is_dir():
        print(f"ERROR: results directory does not exist: {results_dir}")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Results directory : {results_dir}")
    print(f"Output directory  : {output_dir}")
    print(f"Baker systems     : {len(BAKER_SYSTEMS)}")
    print(f"Methods           : {', '.join(METHODS)}")
    print()

    # -- 2D landscape plots --
    print("Generating 2D RMSD landscape plots...")
    landscape_files = generate_landscape_plots(results_dir, output_dir)

    # -- Speedup violin --
    violin_file = None
    if not args.skip_violin:
        print()
        violin_file = run_speedup_violin(output_dir)

    # -- Manifest --
    print_manifest(landscape_files, violin_file)

    if not landscape_files and violin_file is None:
        print("\nNo files were generated. Check that --results-dir is correct.")
        sys.exit(1)


if __name__ == "__main__":
    main()

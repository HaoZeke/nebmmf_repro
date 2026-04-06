#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "rgpycrumbs>=1.4",
#   "chemparseplot>=1.4",
#   "matplotlib>=3.9",
#   "numpy>=2.0",
#   "adjustText",
#   "h5py",
#   "ase>=3.26",
#   "pandas",
#   "seaborn",
# ]
# ///
"""
Regenerate ALL paper figures from current NEB results.

Must be run from the nebmmf root directory with the eongpu pixi env.
Reads from eonRuns/results/03_neb/{system}/{method}/ and writes to
the roneb_tex paper image directories.

Usage:
    pixi run -e eongpu -- python scripts/regen_all_figures.py
"""
import sys
import subprocess
from pathlib import Path

# Paths
NEBMMF = Path(__file__).resolve().parent.parent
EON_RESULTS = NEBMMF / "eonRuns" / "results" / "03_neb"
PAPER_IMGS = Path.home() / "Git/Github/TeX/roneb_tex/v2/FrontiersChemistry/imgs"
SUPPL_IMGS = PAPER_IMGS / "suppl"

BAKER_SYSTEMS = [
    "01_hcn", "02_hcch", "03_h2co", "04_ch3o", "05_cyclopropyl",
    "06_bicyclobutane", "08_formyloxyethyl", "09_parentdielsalder",
    "10_tetrazine", "11_trans_butadiene", "12_ethane_h2_abstraction",
    "13_hf_abstraction", "14_vinyl_alcohol", "15_hocl", "16_h2po4_anion",
    "17_claisen", "18_silylene_insertion", "19_hnccs", "20_hconh3_cation",
    "21_acrolein_rot", "22_hconhoh", "23_hcn_h2", "24_h2cnh", "25_hcnh2",
]

TEAL = "#004D40"
CORAL = "#FF655D"


def run_plt_neb(args: list[str], desc: str) -> bool:
    """Run plt-neb by importing directly (bypasses uvx subprocess issue)."""
    import importlib
    old_argv = sys.argv[:]
    sys.argv = ["plt-neb"] + args
    try:
        mod = importlib.import_module("rgpycrumbs.eon.plt_neb")
        importlib.reload(mod)  # reset state between calls
        mod.main()
        print(f"  OK: {desc}")
        return True
    except Exception as e:
        print(f"  FAIL: {desc}: {e}")
        return False
    finally:
        sys.argv = old_argv


def gen_1d_profile(system: str, method: str, rc_mode: str, output: Path) -> bool:
    """Generate a 1D NEB profile plot."""
    con = EON_RESULTS / system / method / "neb.con"
    dat_pattern = str(EON_RESULTS / system / method / "neb_*.dat")
    if not con.exists():
        return False
    output.parent.mkdir(parents=True, exist_ok=True)
    return run_plt_neb([
        "--con-file", str(con),
        "--output-file", str(output),
        "--plot-type", "profile",
        "--rc-mode", rc_mode,
        "--input-dat-pattern", dat_pattern,
        "--plot-structures", "crit_points",
        "--dpi", "300",
        "--figsize", "5.37", "3.37",
        "--fontsize-base", "12",
        "--zoom-ratio", "0.25",
        "--draw-product", "30.0", "60.0", "0.1",
        "--draw-reactant", "-1.0", "40.0", "0.1",
        "--draw-saddle", "25.0", "40.0", "0.1",
        "--ase-rotation", "0x,90y,0z",
    ], f"{system}/{method} 1D_{rc_mode}")


def gen_2d_landscape(system: str, method: str, output: Path,
                     additional_method: str = None) -> bool:
    """Generate a 2D RMSD landscape plot."""
    con = EON_RESULTS / system / method / "neb.con"
    dat_pattern = str(EON_RESULTS / system / method / "neb_*.dat")
    path_pattern = str(EON_RESULTS / system / method / "neb_path*.con")
    if not con.exists():
        return False
    output.parent.mkdir(parents=True, exist_ok=True)

    args = [
        "--con-file", str(con),
        "--output-file", str(output),
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
    ]

    # Overlay the other method's saddle point if available
    if additional_method:
        sp_con = EON_RESULTS / system / additional_method / "sp.con"
        if sp_con.exists():
            args.extend(["--additional-con", str(sp_con), additional_method.upper()])

    return run_plt_neb(args, f"{system}/{method} 2D_rmsd")


def gen_dumbbell_plot(output: Path) -> bool:
    """Generate the dumbbell/lollipop plot comparing CI-NEB vs OCI-NEB."""
    import re
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # Short chemical labels for each Baker system
    LABELS = {
        "01_hcn": "HCN",
        "02_hcch": "HCCH",
        "03_h2co": "H2CO",
        "04_ch3o": "CH3O",
        "05_cyclopropyl": "Cyclopropyl",
        "06_bicyclobutane": "Bicyclobutane",
        "08_formyloxyethyl": "Formyloxyethyl",
        "09_parentdielsalder": "Diels-Alder",
        "10_tetrazine": "Tetrazine",
        "11_trans_butadiene": "Butadiene",
        "12_ethane_h2_abstraction": "H2 abstraction",
        "13_hf_abstraction": "HF abstraction",
        "14_vinyl_alcohol": "Vinyl alcohol",
        "15_hocl": "HOCl",
        "16_h2po4_anion": "H2PO4-",
        "17_claisen": "Claisen",
        "18_silylene_insertion": "Silylene ins.",
        "19_hnccs": "HNCCS",
        "20_hconh3_cation": "HCONH3+",
        "21_acrolein_rot": "Acrolein rot.",
        "22_hconhoh": "HCONHOH",
        "23_hcn_h2": "HCN+H2",
        "24_h2cnh": "H2CNH",
        "25_hcnh2": "HCNH2",
    }

    systems = []
    ci_calls = []
    mmf_calls = []

    for sys in BAKER_SYSTEMS:
        ci_file = EON_RESULTS / sys / "cineb" / "results.dat"
        mmf_file = EON_RESULTS / sys / "mmf" / "results.dat"
        if not ci_file.exists() or not mmf_file.exists():
            continue
        ci = int(re.search(r"(\d+) total_force_calls", ci_file.read_text()).group(1))
        mmf = int(re.search(r"(\d+) total_force_calls", mmf_file.read_text()).group(1))
        label = LABELS.get(sys, sys.split("_", 1)[1])
        systems.append(label)
        ci_calls.append(ci)
        mmf_calls.append(mmf)

    # Sort by CI-NEB cost (ascending)
    order = sorted(range(len(ci_calls)), key=lambda i: ci_calls[i])
    systems = [systems[i] for i in order]
    ci_calls = [ci_calls[i] for i in order]
    mmf_calls = [mmf_calls[i] for i in order]

    fig, ax = plt.subplots(figsize=(8, 7))
    y = np.arange(len(systems))

    # Dumbbell lines
    for i in range(len(systems)):
        ax.plot([mmf_calls[i], ci_calls[i]], [y[i], y[i]], color="#BBBBBB", lw=1.5, alpha=0.7)

    ax.scatter(ci_calls, y, color=TEAL, s=50, zorder=5, label="CI-NEB")
    ax.scatter(mmf_calls, y, color=CORAL, s=50, zorder=5, label="OCI-NEB")

    ax.set_yticks(y)
    ax.set_yticklabels(systems, fontsize=10, color=TEAL)
    ax.set_xlabel("Force evaluations", fontsize=13, color=TEAL)
    ax.legend(fontsize=11, loc="lower right")
    ax.invert_yaxis()

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_color(TEAL)
    ax.tick_params(colors=TEAL)
    ax.grid(axis="x", alpha=0.15)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(str(output), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  OK: dumbbell_plot")
    return True


def gen_speedup_violin(output: Path) -> bool:
    """Generate the speedup ratio violin plot."""
    import re
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    ratios = []
    for sys in BAKER_SYSTEMS:
        ci_file = EON_RESULTS / sys / "cineb" / "results.dat"
        mmf_file = EON_RESULTS / sys / "mmf" / "results.dat"
        if not ci_file.exists() or not mmf_file.exists():
            continue
        ci = int(re.search(r"(\d+) total_force_calls", ci_file.read_text()).group(1))
        mmf = int(re.search(r"(\d+) total_force_calls", mmf_file.read_text()).group(1))
        ratios.append(ci / mmf)

    fig, ax = plt.subplots(figsize=(3.5, 5))

    parts = ax.violinplot(ratios, positions=[0], showmeans=False,
                          showmedians=False, showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor(CORAL)
        pc.set_edgecolor(TEAL)
        pc.set_alpha(0.6)
        pc.set_linewidth(1.2)

    bp = ax.boxplot(ratios, positions=[0], widths=0.15, patch_artist=True,
                    boxprops=dict(facecolor="white", edgecolor=TEAL, linewidth=1.2),
                    whiskerprops=dict(color=TEAL, linewidth=1.2),
                    capprops=dict(color=TEAL, linewidth=1.2),
                    medianprops=dict(color=CORAL, linewidth=2),
                    flierprops=dict(marker="o", markerfacecolor=TEAL,
                                    markeredgecolor=TEAL, markersize=4))

    jitter = np.random.default_rng(42).uniform(-0.06, 0.06, size=len(ratios))
    ax.scatter(jitter, ratios, color=TEAL, s=18, alpha=0.7, zorder=5, edgecolors="none")
    ax.axhline(y=1.0, color="#999999", linestyle="--", linewidth=0.8, zorder=1)

    ax.set_ylabel("Speedup Ratio (CI-NEB / OCI-NEB)", fontsize=11, color=TEAL)
    ax.set_xticks([0])
    ax.set_xticklabels([f"Baker Benchmark\n({len(ratios)} systems)"],
                       fontsize=10, color=TEAL)
    ax.tick_params(axis="y", colors=TEAL, labelsize=10)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["bottom", "left"]:
        ax.spines[spine].set_color(TEAL)
        ax.spines[spine].set_linewidth(1.2)

    ax.set_xlim(-0.6, 0.6)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(str(output), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    import statistics
    print(f"  OK: speedup_violin (median={statistics.median(ratios):.2f}x, "
          f"mean={statistics.mean(ratios):.2f}x)")
    return True


def main():
    print("=" * 60)
    print("  Regenerating ALL paper figures")
    print("=" * 60)
    generated = []
    failed = []

    # --- Main paper figures ---
    print("\n[1/6] Dumbbell plot (Fig 1)")
    if gen_dumbbell_plot(PAPER_IMGS / "dumbbell_plot.png"):
        generated.append("dumbbell_plot.png")
    else:
        failed.append("dumbbell_plot.png")

    print("\n[2/6] Case study: HNCCS (System 19)")
    for method, suffix in [("cineb", "cineb"), ("mmf", "mmf")]:
        for rc in ["path"]:
            out = PAPER_IMGS / f"hnccs_{suffix}_1D_{rc}.png"
            if gen_1d_profile("19_hnccs", method, rc, out):
                generated.append(out.name)
            else:
                failed.append(out.name)

    out = PAPER_IMGS / "hnccs_mmf_2D.png"
    if gen_2d_landscape("19_hnccs", "mmf", out, additional_method="cineb"):
        generated.append(out.name)
    else:
        failed.append(out.name)

    print("\n[3/6] Case study: HCONH3+ (System 20)")
    out = PAPER_IMGS / "cases" / "hconh3cat_cineb_2D.png"
    if gen_2d_landscape("20_hconh3_cation", "cineb", out, additional_method="mmf"):
        generated.append(out.name)
    else:
        failed.append(out.name)

    print("\n[4/6] Case study: Claisen (System 17)")
    out = PAPER_IMGS / "cases" / "17_claisen_static_neb_2D.png"
    if gen_2d_landscape("17_claisen", "mmf", out, additional_method="cineb"):
        generated.append(out.name)
    else:
        failed.append(out.name)

    # --- Supplementary figures ---
    print("\n[5/6] Speedup violin plot")
    if gen_speedup_violin(SUPPL_IMGS / "speedup_ratio_violin.png"):
        generated.append("speedup_ratio_violin.png")
    else:
        failed.append("speedup_ratio_violin.png")

    print("\n[6/6] Baker 2D landscapes (all 24 systems)")
    baker_2d_dir = SUPPL_IMGS / "baker_2d"
    for sys in BAKER_SYSTEMS:
        for method in ["cineb", "mmf"]:
            out = baker_2d_dir / sys / f"{method}_landscape.png"
            other = "mmf" if method == "cineb" else "cineb"
            if gen_2d_landscape(sys, method, out, additional_method=other):
                generated.append(f"baker_2d/{sys}/{method}")
            else:
                failed.append(f"baker_2d/{sys}/{method}")

    # --- Summary ---
    print("\n" + "=" * 60)
    print(f"  Generated: {len(generated)}")
    print(f"  Failed:    {len(failed)}")
    if failed:
        print("  Failures:")
        for f in failed:
            print(f"    - {f}")
    print("=" * 60)


if __name__ == "__main__":
    main()

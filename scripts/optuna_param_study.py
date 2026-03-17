#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "optuna",
#   "plotly",
#   "kaleido",
#   "rgpycrumbs @ git+https://github.com/HaoZeke/rgpycrumbs@eonMLFlow",
# ]
# ///

"""
Optuna-based parameter sensitivity study for OCI-NEB (MMF) parameters.

Runs a TPE sampler over the 7 OCI-NEB hyperparameters on 4 representative
Baker systems, minimizing total force evaluations. After the study, fANOVA
importance analysis identifies which parameters actually drive performance.
"""

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import optuna
from optuna.importance import get_param_importances

from rgpycrumbs.eon.helpers import write_eon_config

CALLS_REGEX = re.compile(r"destroyed after (\d+) calls")

REPRESENTATIVE_SYSTEMS = [
    "05_cyclopropyl",
    "17_claisen",
    "19_hnccs",
    "20_hconh3_cation",
]

# Large penalty for failed/non-converged runs
PENALTY_FORCE_CALLS = 50000


def find_model_path(base_dir: Path) -> Path:
    """Auto-detect the .pt model file under results/00_models/."""
    models_dir = base_dir / "results" / "00_models"
    pts = list(models_dir.glob("*.pt"))
    if not pts:
        raise FileNotFoundError(f"No .pt model found in {models_dir}")
    # Pick the first (there should be exactly one)
    return pts[0]


def parse_force_calls(run_dir: Path) -> int | None:
    """Parse total force calls from eonclient log files."""
    for log in run_dir.glob("*.log"):
        try:
            txt = log.read_text(errors="ignore")
            m = CALLS_REGEX.search(txt)
            if m:
                return int(m.group(1))
        except Exception:
            pass
    return None


def check_convergence(run_dir: Path) -> bool:
    """Check whether the NEB run converged (status 0 in results.dat)."""
    results_file = run_dir / "results.dat"
    if not results_file.exists():
        return False
    txt = results_file.read_text(errors="ignore")
    first_line = txt.splitlines()[0] if txt.strip() else ""
    # termination_reason line starts with the status code
    if "termination_reason" in first_line:
        code = first_line.split()[0]
        try:
            return int(code) == 0
        except (ValueError, TypeError):
            return False
    return False


def run_single_system(
    system: str,
    params: dict,
    base_dir: Path,
    model_path: Path,
    results_dir: str,
    n_images: int = 8,
) -> int:
    """Run eonclient for one system with the given OCI-NEB params.

    Returns the number of force calls, or PENALTY_FORCE_CALLS on failure.
    """
    endpoints_dir = base_dir / "results" / "01_endpoints" / system
    idpp_dir = base_dir / "results" / "02_idpp" / system

    reactant = endpoints_dir / "reactant.con"
    product = endpoints_dir / "product.con"
    idpp_path = idpp_dir / "idppPath.dat"

    for f in [reactant, product, idpp_path]:
        if not f.exists():
            print(f"  [WARN] Missing input: {f}")
            return PENALTY_FORCE_CALLS

    # Build NEB parameters matching the snakemake do_neb rule
    neb_params = {
        "images": n_images,
        "energy_weighted": "true",
        "ew_ksp_min": 0.972,
        "ew_ksp_max": 9.72,
        "ew_trigger": 0.5,
        "initializer": "sidpp",
        "sidpp_growth_alpha": 0.33,
        "oversampling": "false",
        "oversampling_factor": 8,
        "minimize_endpoints": "false",
        "climbing_image_method": "true",
        "climbing_image_converged_only": "true",
        "ci_after": 0.5,
        "ci_after_rel": 0.8,
        # OCI-NEB (MMF) parameters from the trial
        "ci_mmf": "true",
        "ci_mmf_after": params["ci_mmf_after"],
        "ci_mmf_after_rel": params["ci_mmf_after_rel"],
        "ci_mmf_angle": params["ci_mmf_angle"],
        "ci_mmf_penalty_strength": params["ci_mmf_penalty_strength"],
        "ci_mmf_penalty_base": params["ci_mmf_penalty_base"],
        "ci_mmf_nsteps": params["ci_mmf_nsteps"],
        "ci_mmf_ci_stability_count": params["ci_mmf_ci_stability_count"],
    }

    neb_settings = {
        "Main": {"job": "nudged_elastic_band", "random_seed": 706253457},
        "Potential": {"potential": "metatomic"},
        "Metatomic": {
            "model_path": str(model_path.absolute()),
            "device": "cuda",
        },
        "Nudged Elastic Band": neb_params,
        "Dimer": {
            "improved": "true",
            "opt_method": "cg",
            "remove_rotation": "false",
            "converged_angle": 10.0,
        },
        "Optimizer": {
            "max_iterations": 1000,
            "opt_method": "lbfgs",
            "max_move": 0.1,
            "converged_force": 0.05,
        },
    }

    with tempfile.TemporaryDirectory(prefix=f"optuna_{system}_") as tmp:
        work_dir = Path(tmp)
        write_eon_config(work_dir, neb_settings)
        shutil.copy2(reactant, work_dir / "reactant.con")
        shutil.copy2(product, work_dir / "product.con")
        shutil.copy2(idpp_path, work_dir / "idppPath.dat")

        # Copy IDPP path images if they exist
        idpp_path_dir = idpp_dir / "path"
        if idpp_path_dir.is_dir():
            dest_path_dir = work_dir / "path"
            dest_path_dir.mkdir(exist_ok=True)
            for img_file in sorted(idpp_path_dir.glob("*.con")):
                shutil.copy2(img_file, dest_path_dir / img_file.name)

        try:
            subprocess.run(
                ["eonclient"],
                cwd=work_dir,
                check=True,
                capture_output=True,
                timeout=600,
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
            print(f"  [WARN] eonclient failed for {system}: {exc}")
            return PENALTY_FORCE_CALLS

        if not check_convergence(work_dir):
            # Non-converged: still try to get force calls, but add penalty
            calls = parse_force_calls(work_dir)
            if calls is not None:
                return calls + PENALTY_FORCE_CALLS // 2
            return PENALTY_FORCE_CALLS

        calls = parse_force_calls(work_dir)
        if calls is None:
            return PENALTY_FORCE_CALLS
        return calls


def objective(
    trial: optuna.Trial,
    base_dir: Path,
    model_path: Path,
    results_dir: str,
) -> float:
    """Optuna objective: sum of force calls across 4 representative systems."""
    params = {
        "ci_mmf_after_rel": trial.suggest_float(
            "ci_mmf_after_rel", 0.2, 0.8
        ),
        "ci_mmf_angle": trial.suggest_float("ci_mmf_angle", 0.3, 1.5),
        "ci_mmf_penalty_strength": trial.suggest_float(
            "ci_mmf_penalty_strength", 0.1, 5.0
        ),
        "ci_mmf_penalty_base": trial.suggest_float(
            "ci_mmf_penalty_base", 0.1, 1.0
        ),
        "ci_mmf_after": trial.suggest_float("ci_mmf_after", 0.05, 0.5),
        "ci_mmf_nsteps": trial.suggest_int("ci_mmf_nsteps", 100, 2000),
        "ci_mmf_ci_stability_count": trial.suggest_int(
            "ci_mmf_ci_stability_count", 2, 10
        ),
    }

    total_calls = 0
    for system in REPRESENTATIVE_SYSTEMS:
        calls = run_single_system(
            system, params, base_dir, model_path, results_dir
        )
        total_calls += calls
        trial.report(total_calls, step=REPRESENTATIVE_SYSTEMS.index(system))
        if trial.should_prune():
            raise optuna.TrialPruned()

    return total_calls


def main():
    parser = argparse.ArgumentParser(
        description="Optuna parameter sensitivity study for OCI-NEB"
    )
    parser.add_argument(
        "--n-trials",
        type=int,
        default=100,
        help="Number of optuna trials (default: 100)",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Parallel jobs (default: 1 for sequential GPU runs)",
    )
    parser.add_argument(
        "--results-dir",
        type=str,
        default="eonRuns/results/03_neb",
        help="NEB results directory (default: eonRuns/results/03_neb)",
    )
    parser.add_argument(
        "--model-path",
        type=str,
        default=None,
        help="Path to .pt model file (default: auto-detect from config)",
    )
    args = parser.parse_args()

    # Resolve base directory (the eonRuns parent)
    script_dir = Path(__file__).resolve().parent
    base_dir = script_dir.parent / "eonRuns"
    if not base_dir.is_dir():
        base_dir = Path.cwd()
        if (base_dir / "eonRuns").is_dir():
            base_dir = base_dir / "eonRuns"

    if args.model_path:
        model_path = Path(args.model_path).resolve()
    else:
        model_path = find_model_path(base_dir)
    print(f"Using model: {model_path}")

    # Set up results output directory
    out_dir = base_dir.parent / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    db_path = out_dir / "optuna_study.db"
    storage = f"sqlite:///{db_path}"

    study = optuna.create_study(
        study_name="oci_neb_param_sensitivity",
        storage=storage,
        direction="minimize",
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(seed=42),
        pruner=optuna.pruners.MedianPruner(
            n_startup_trials=10, n_warmup_steps=1
        ),
    )

    study.optimize(
        lambda trial: objective(trial, base_dir, model_path, args.results_dir),
        n_trials=args.n_trials,
        n_jobs=args.n_jobs,
        show_progress_bar=True,
    )

    # -- Parameter importance (fANOVA) --
    print("\n" + "=" * 60)
    print("Parameter Importances (fANOVA)")
    print("=" * 60)
    importances = get_param_importances(study)
    for param, importance in importances.items():
        bar = "#" * int(importance * 40)
        print(f"  {param:30s}  {importance:.4f}  {bar}")

    # -- Best parameters --
    print("\n" + "=" * 60)
    print("Best Parameters")
    print("=" * 60)
    print(f"  Total force calls: {study.best_value}")
    for k, v in study.best_params.items():
        print(f"  {k:30s}  {v}")

    # -- Visualization --
    try:
        fig = optuna.visualization.plot_param_importances(study)
        plot_path = out_dir / "optuna_importance.png"
        fig.write_image(str(plot_path), width=900, height=500, scale=2)
        print(f"\nImportance plot saved to: {plot_path}")
    except Exception as exc:
        print(f"\n[WARN] Could not save importance plot: {exc}")

    print(f"Study database saved to: {db_path}")


if __name__ == "__main__":
    main()

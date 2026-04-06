#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "optuna",
#   "matplotlib",
#   "scikit-learn",
#   "pyyaml",
#   "rgpycrumbs",
# ]
# ///

"""
Optuna 2-parameter study for OCI-NEB: after_rel + angle only.

Fixes nsteps=1000 and stability_count=3 (insensitive per 5-param fANOVA),
isolating the two parameters that matter most.  The 2D contour of this
study shows angle dominates at 81% when the other three are held constant.

Usage:
    CUBLAS_WORKSPACE_CONFIG=:4096:8 pixi run -e eongpu \
        python scripts/optuna_2param.py --n-trials 200
"""

import argparse
import re
import shutil
import subprocess
import tempfile
import threading
from pathlib import Path

import optuna
import yaml
from optuna.importance import get_param_importances

from rgpycrumbs.eon.helpers import write_eon_config

CALLS_REGEX = re.compile(r"(\d+) total_force_calls")

BAKER_SYSTEMS = [
    "13_hf_abstraction",
    "09_parentdielsalder",
    "18_silylene_insertion",
    "12_ethane_h2_abstraction",
    "17_claisen",
    "06_bicyclobutane",
    "11_trans_butadiene",
    "08_formyloxyethyl",
    "01_hcn",
    "02_hcch",
    "03_h2co",
    "04_ch3o",
    "05_cyclopropyl",
    "10_tetrazine",
    "14_vinyl_alcohol",
    "15_hocl",
    "16_h2po4_anion",
    "19_hnccs",
    "20_hconh3_cation",
    "21_acrolein_rot",
    "22_hconhoh",
    "23_hcn_h2",
    "24_h2cnh",
    "25_hcnh2",
]

PENALTY_FORCE_CALLS = 50000

# Fixed parameters (insensitive per 5-param fANOVA)
FIXED_NSTEPS = 1000
FIXED_AFTER = 0.3
FIXED_STABILITY_COUNT = 3


def find_model_path(base_dir: Path) -> Path:
    models_dir = base_dir / "results" / "00_models"
    pts = list(models_dir.glob("*.pt"))
    if not pts:
        raise FileNotFoundError(f"Found no .pt model in {models_dir}")
    return pts[0]


def parse_eon_results(run_dir: Path) -> tuple[bool, int | None]:
    results_file = run_dir / "results.dat"
    if not results_file.exists():
        return False, None

    converged = False
    calls = None
    try:
        txt = results_file.read_text(errors="ignore")
        lines = txt.splitlines()

        if lines and "termination_reason" in lines[0]:
            code = lines[0].split()[0]
            converged = (code == "0")

        m = CALLS_REGEX.search(txt)
        if m:
            calls = int(m.group(1))
    except (OSError, ValueError):
        pass

    return converged, calls


def run_single_system(
    system: str,
    params: dict,
    base_dir: Path,
    model_path: Path,
    n_images: int = 8,
) -> int:
    endpoints_dir = base_dir / "results" / "01_endpoints" / system
    reactant = endpoints_dir / "reactant.con"
    product = endpoints_dir / "product.con"

    for f in [reactant, product]:
        if not f.exists():
            print(f"  [WARN] Missing input structure: {f}")
            return PENALTY_FORCE_CALLS

    neb_params = {
        "images": n_images,
        "energy_weighted": "true",
        "ew_ksp_min": 0.972,
        "ew_ksp_max": 9.72,
        "ew_trigger": 0.5,
        "initializer": "sidpp",
        "oversampling": "false",
        "oversampling_factor": 8,
        "minimize_endpoints": "false",
        "climbing_image_method": "true",
        "climbing_image_converged_only": "true",
        "ci_after": 0.5,
        "ci_after_rel": 0.8,
        "ci_mmf": "true",
        "ci_mmf_after": FIXED_AFTER,
        "ci_mmf_after_rel": params["ci_mmf_after_rel"],
        "ci_mmf_angle": params["ci_mmf_angle"],
        "ci_mmf_nsteps": FIXED_NSTEPS,
        "ci_mmf_ci_stability_count": FIXED_STABILITY_COUNT,
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

    ramdisk = Path("/dev/shm")
    base_tmp = ramdisk if ramdisk.exists() and ramdisk.is_dir() else None

    with tempfile.TemporaryDirectory(prefix=f"opt2p_{system}_", dir=base_tmp) as tmp:
        work_dir = Path(tmp)
        write_eon_config(work_dir, neb_settings)
        shutil.copy2(reactant, work_dir / "reactant.con")
        shutil.copy2(product, work_dir / "product.con")

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

        converged, calls = parse_eon_results(work_dir)

        if not converged:
            if calls is not None:
                return calls + PENALTY_FORCE_CALLS // 2
            return PENALTY_FORCE_CALLS

        if calls is None:
            return PENALTY_FORCE_CALLS
        return calls


_CI_BASELINES: dict[str, int] = {}
_CI_LOCK = threading.Lock()


def run_ci_neb(system: str, base_dir: Path, model_path: Path, n_images: int = 8) -> int:
    endpoints_dir = base_dir / "results" / "01_endpoints" / system
    reactant = endpoints_dir / "reactant.con"
    product = endpoints_dir / "product.con"

    for f in [reactant, product]:
        if not f.exists():
            return PENALTY_FORCE_CALLS

    neb_params = {
        "images": n_images,
        "energy_weighted": "true",
        "ew_ksp_min": 0.972,
        "ew_ksp_max": 9.72,
        "ew_trigger": 0.5,
        "initializer": "sidpp",
        "oversampling": "false",
        "oversampling_factor": 8,
        "minimize_endpoints": "false",
        "climbing_image_method": "true",
        "climbing_image_converged_only": "true",
        "ci_after": 0.5,
        "ci_after_rel": 0.8,
        "ci_mmf": "false",
    }

    neb_settings = {
        "Main": {"job": "nudged_elastic_band", "random_seed": 706253457},
        "Potential": {"potential": "metatomic"},
        "Metatomic": {
            "model_path": str(model_path.absolute()),
            "device": "cuda",
        },
        "Nudged Elastic Band": neb_params,
        "Optimizer": {
            "max_iterations": 1000,
            "opt_method": "lbfgs",
            "max_move": 0.1,
            "converged_force": 0.05,
        },
    }

    ramdisk = Path("/dev/shm")
    base_tmp = ramdisk if ramdisk.exists() and ramdisk.is_dir() else None

    with tempfile.TemporaryDirectory(prefix=f"cineb_{system}_", dir=base_tmp) as tmp:
        work_dir = Path(tmp)
        write_eon_config(work_dir, neb_settings)
        shutil.copy2(reactant, work_dir / "reactant.con")
        shutil.copy2(product, work_dir / "product.con")

        try:
            subprocess.run(
                ["eonclient"],
                cwd=work_dir,
                check=True,
                capture_output=True,
                timeout=600,
            )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            return PENALTY_FORCE_CALLS

        converged, calls = parse_eon_results(work_dir)
        return calls if calls is not None else PENALTY_FORCE_CALLS


def get_ci_baseline(system: str, base_dir: Path, model_path: Path) -> int:
    with _CI_LOCK:
        if system in _CI_BASELINES:
            return _CI_BASELINES[system]

    calls = run_ci_neb(system, base_dir, model_path)

    with _CI_LOCK:
        if system not in _CI_BASELINES:
            _CI_BASELINES[system] = calls
            print(f"  [CI baseline] {system}: {calls}")
        return _CI_BASELINES[system]


def objective(
    trial: optuna.Trial,
    base_dir: Path,
    model_path: Path,
) -> float:
    """2-parameter objective: after_rel and angle only."""
    params = {
        "ci_mmf_after_rel": trial.suggest_float("ci_mmf_after_rel", 0.2, 0.8),
        "ci_mmf_angle": trial.suggest_float("ci_mmf_angle", 0.5, 1.0),
    }

    print(f"\n[Trial {trial.number}] Parameters: {params}")

    total_calls = 0
    regression_penalty = 0
    for idx, system in enumerate(BAKER_SYSTEMS):
        calls = run_single_system(system, params, base_dir, model_path)
        ci_baseline = get_ci_baseline(system, base_dir, model_path)

        ratio = calls / max(ci_baseline, 1)
        if ratio > 1.0:
            regression_penalty += (calls - ci_baseline) * 100
            print(
                f"  [Trial {trial.number} | REGRESSION] {system}: ratio={ratio:.2f} "
                f"(OCI-NEB={calls}, CI-NEB={ci_baseline})"
            )

        total_calls += calls
        trial.report(total_calls + regression_penalty, step=idx)
        if trial.should_prune():
            print(f"  [Trial {trial.number}] Pruned at {system}.")
            raise optuna.TrialPruned()

    print(f"[Trial {trial.number}] Completed: {total_calls + regression_penalty}")
    return total_calls + regression_penalty


def main():
    parser = argparse.ArgumentParser(
        description="Optuna 2-parameter study for OCI-NEB (after_rel + angle)"
    )
    parser.add_argument("--n-trials", type=int, default=200)
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    base_dir = script_dir.parent / "eonRuns"
    if not base_dir.is_dir():
        base_dir = Path.cwd()
        if (base_dir / "eonRuns").is_dir():
            base_dir = base_dir / "eonRuns"

    model_path = find_model_path(base_dir)
    print(f"Model: {model_path}")

    out_dir = base_dir.parent / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    db_path = out_dir / "optuna_2param.db"
    storage = f"sqlite:///{db_path}"

    study = optuna.create_study(
        study_name="oci_neb_2param",
        storage=storage,
        direction="minimize",
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(seed=42, n_startup_trials=20),
        pruner=optuna.pruners.MedianPruner(n_startup_trials=25, n_warmup_steps=1),
    )

    study.optimize(
        lambda trial: objective(trial, base_dir, model_path),
        n_trials=args.n_trials,
        n_jobs=args.n_jobs,
        show_progress_bar=True,
    )

    print("\n" + "=" * 60)
    print("Parameter Importances (fANOVA)")
    print("=" * 60)
    importances = get_param_importances(study)
    for param, importance in importances.items():
        bar = "#" * int(importance * 40)
        print(f"  {param:30s}  {importance:.4f}  {bar}")

    print("\n" + "=" * 60)
    print("Best Parameters")
    print("=" * 60)
    print(f"  Total force evaluations: {study.best_value}")
    for k, v in study.best_params.items():
        print(f"  {k:30s}  {v}")

    try:
        import matplotlib.pyplot as plt
        from optuna.visualization.matplotlib import plot_param_importances, plot_contour

        plt.rcParams.update({"font.size": 18, "axes.labelsize": 20,
                             "axes.titlesize": 22, "xtick.labelsize": 16,
                             "ytick.labelsize": 16})

        fig_imp = plot_param_importances(study)
        fig_imp.set_size_inches(10, 6)
        plt.tight_layout()
        imp_path = out_dir / "optuna_2param_importance.png"
        plt.savefig(imp_path, dpi=300, bbox_inches="tight")
        print(f"\nSaved importance plot: {imp_path}")
        plt.close()

        fig_contour = plot_contour(study, params=["ci_mmf_after_rel", "ci_mmf_angle"])
        fig_contour.set_size_inches(10, 8)
        plt.tight_layout()
        contour_path = out_dir / "optuna_2param_contour.png"
        plt.savefig(contour_path, dpi=300, bbox_inches="tight")
        print(f"Saved contour plot: {contour_path}")
        plt.close()

    except Exception as exc:
        print(f"\n[WARN] Failed writing plots: {exc}")

    print(f"Study database: {db_path}")

    best = study.best_params
    params_yaml = out_dir / "optuna_2param_best.yaml"
    best_config = {
        "ci_mmf_after_rel": round(best["ci_mmf_after_rel"], 4),
        "ci_mmf_angle": round(best["ci_mmf_angle"], 4),
        "fixed": {
            "ci_mmf_nsteps": FIXED_NSTEPS,
            "ci_mmf_after": FIXED_AFTER,
            "ci_mmf_ci_stability_count": FIXED_STABILITY_COUNT,
        },
    }
    with open(params_yaml, "w") as f:
        yaml.dump(best_config, f, default_flow_style=False)
    print(f"Exported best params: {params_yaml}")


if __name__ == "__main__":
    main()

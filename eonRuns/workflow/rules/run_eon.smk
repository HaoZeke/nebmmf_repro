# -*- mode:snakemake; -*-
"""
NEB optimization using eOn with CI-NEB and optional OCI-NEB (MMF) refinement.

Uses rgpycrumbs.eon.helpers.write_eon_config to generate eOn configuration
programmatically from YAML config, following the atomistic cookbook pattern.
"""

from pathlib import Path
from rgpycrumbs.eon.helpers import write_eon_config
import shutil
import subprocess
import os

_pet_type = config.get("pet_mad", {}).get("type", "pet-mad")
_pet_version = config["pet_mad"]["version"]
_pet_name = f"{_pet_type}-{_pet_version}"

# NEB method settings: cineb (baseline) vs mmf (OCI-NEB)
_METHOD_OVERRIDES = {
    "cineb": {
        "ci_mmf": "false",
    },
    "mmf": {
        "ci_mmf": "true",
        "ci_mmf_after": 0.1,
        "ci_mmf_after_rel": 0.5,
        "ci_mmf_penalty_strength": 1.5,
        "ci_mmf_penalty_base": 0.4,
        "ci_mmf_angle": 0.9,
        "ci_mmf_nsteps": 1000,
    },
    "strict": {
        "ci_mmf": "true",
        "ci_mmf_after": 0.1,
        "ci_mmf_after_rel": 0.7,
        "ci_mmf_penalty_strength": 0.5,
        "ci_mmf_penalty_base": 0.4,
        "ci_mmf_angle": 0.5,
        "ci_mmf_nsteps": 1000,
    },
}

# Parameter sensitivity analysis is handled by scripts/optuna_param_study.py
# using optuna TPE + fANOVA importance, replacing the earlier grid search approach.


rule do_minimization:
    input:
        endpoint=f"{config['paths']['endpoints']}/{{system}}/{{endpoint}}_pre_aligned.con",
        model=f"{config['paths']['models']}/{_pet_name}.pt",
    output:
        endpoint=f"{config['paths']['endpoints']}/{{system}}/{{endpoint}}_minimized.con",
    run:
        min_settings = {
            "Main": {"job": "minimization", "random_seed": 706253457},
            "Potential": {"potential": "metatomic"},
            "Metatomic": {
                "model_path": str(Path(input.model).absolute()),
                "device": "cuda",
            },
            "Optimizer": {
                "max_iterations": 2000,
                "opt_method": "lbfgs",
                "max_move": 0.1,
                "converged_force": 0.01,
            },
        }

        work_dir = Path(output.endpoint).parent / f"_min_{wildcards.endpoint}"
        work_dir.mkdir(parents=True, exist_ok=True)

        write_eon_config(work_dir, min_settings)
        shutil.copy2(os.path.abspath(input.endpoint), work_dir / "pos.con")

        subprocess.run(["eonclient"], cwd=work_dir, check=True)
        shutil.copy2(work_dir / "min.con", output.endpoint)


rule do_neb:
    input:
        idpp_path=f"{config['paths']['idpp']}/{{system}}/idppPath.dat",
        reactant=f"{config['paths']['endpoints']}/{{system}}/reactant.con",
        product=f"{config['paths']['endpoints']}/{{system}}/product.con",
        path_images=lambda w: expand(
            f"{config['paths']['idpp']}/{w.system}/path/{{num:02d}}.con",
            num=range(config["common"]["number_of_intermediate_imgs"] + 2),
        ),
        model=f"{config['paths']['models']}/{_pet_name}.pt",
    output:
        results_dat=f"{config['paths']['neb']}/{{system}}/{{method}}/results.dat",
        neb_con=f"{config['paths']['neb']}/{{system}}/{{method}}/neb.con",
        neb_dat=f"{config['paths']['neb']}/{{system}}/{{method}}/neb.dat",
    params:
        opath=f"{config['paths']['neb']}/{{system}}/{{method}}",
    run:
        method = wildcards.method
        overrides = _METHOD_OVERRIDES.get(method, _METHOD_OVERRIDES["cineb"])

        neb_params = {
            "images": config["common"]["number_of_intermediate_imgs"],
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
        }
        neb_params.update(overrides)

        neb_settings = {
            "Main": {"job": "nudged_elastic_band", "random_seed": 706253457},
            "Potential": {"potential": "metatomic"},
            "Metatomic": {
                "model_path": str(Path(input.model).absolute()),
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
            "Debug": {"write_movies": "true"},
        }

        out_path = Path(params.opath)
        out_path.mkdir(parents=True, exist_ok=True)

        write_eon_config(out_path, neb_settings)
        shutil.copy2(os.path.abspath(input.reactant), out_path / "reactant.con")
        shutil.copy2(os.path.abspath(input.product), out_path / "product.con")
        shutil.copy2(os.path.abspath(input.idpp_path), out_path / "idppPath.dat")

        subprocess.run(["eonclient"], cwd=out_path, check=True, capture_output=True)

# -*- mode:snakemake; -*-
"""
PET-MAD model retrieval and conversion.

Downloads and exports model checkpoints from HuggingFace using mtt export.
Supports pet-mad, pet-mad-xs, pet-mad-s variants via config['pet_mad']['type'].
"""

_pet_type = config.get("pet_mad", {}).get("type", "pet-mad")
_pet_version = config["pet_mad"]["version"]
_pet_name = f"{_pet_type}-{_pet_version}"


rule download_and_export_model:
    output:
        protected(f"{config['paths']['models']}/{_pet_name}.pt"),
    params:
        repo="lab-cosmo/upet",
        ckpt_path=f"models/{_pet_name}.ckpt",
    shell:
        """
        mtt export {params.repo} {params.ckpt_path} -o {output}
        """

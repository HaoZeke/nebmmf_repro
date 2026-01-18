# -*- mode:snakemake; -*-

import ase.io
from rgpycrumbs.geom.api.alignment import (
    align_structure_robust,
    IRAConfig,
    AlignmentMethod,
)
import numpy as np
import copy


def align_structures(ref_atm, target_atm, use_ira=True, kmax=1.8):
    """
    Core logic to align target_atm to ref_atm using the robust API.
    """
    config = IRAConfig(enabled=use_ira, kmax=kmax)
    result = align_structure_robust(ref_atm, target_atm, config)

    if result.method == AlignmentMethod.IRA_PERMUTATION:
        print(f"IRA matching successful for system.")
    else:
        print("Standard ASE rotation/translation minimization complete (Fallback).")

    return result.atoms


rule prepare_endpoints_pre:
    """
    Normalizes and aligns initial reactant and product geometries.
    """
    input:
        reactant=lambda wildcards: config["systems"][wildcards.system]["reactant"],
        product=lambda wildcards: config["systems"][wildcards.system]["product"],
    output:
        reactant=f"{config['paths']['endpoints']}/{{system}}/reactant_pre_aligned.con",
        product=f"{config['paths']['endpoints']}/{{system}}/product_pre_aligned.con",
    run:
        reactant_atm = ase.io.read(input.reactant)
        product_atm = ase.io.read(input.product)

        # Standardize cell size and center coordinates for the potential evaluator
        for atm in [reactant_atm, product_atm]:
            atm.set_cell([25, 25, 25])
            atm.center()

        sys_cfg = config["systems"][wildcards.system]
        use_ira = sys_cfg.get("use_ira", True)

        print(f"Processing PRE alignment for system: {wildcards.system}")
        final_product = align_structures(reactant_atm, product_atm, use_ira=use_ira)

        # Calculate initial RMSD for the log
        diff_sq = (reactant_atm.get_positions() - final_product.get_positions()) ** 2
        rmsd = np.sqrt(np.mean(np.sum(diff_sq, axis=1)))
        print(f"Pre-minimization RMSD: {rmsd:.6f} Å")

        ase.io.write(output.reactant, reactant_atm)
        ase.io.write(output.product, final_product)


rule prepare_endpoints_post:
    """
    Aligns minimized endpoints to ensure consistent mapping before NEB.
    """
    input:
        reactant=f"{config['paths']['endpoints']}/{{system}}/reactant_minimized.con",
        product=f"{config['paths']['endpoints']}/{{system}}/product_minimized.con",
    output:
        reactant=f"{config['paths']['endpoints']}/{{system}}/reactant.con",
        product=f"{config['paths']['endpoints']}/{{system}}/product.con",
    run:
        try:
            reactant_atm = ase.io.read(input.reactant)
            product_atm = ase.io.read(input.product)
        except FileNotFoundError as e:
            raise Exception(f"Error: {e}. Missing minimized endpoints.")

        sys_cfg = config["systems"][wildcards.system]
        use_ira = sys_cfg.get("use_ira", True)

        print(f"Processing POST alignment for system: {wildcards.system}")
        final_product = align_structures(reactant_atm, product_atm, use_ira=use_ira)

        # Calculate final RMSD for verification of the mapping
        diff_sq = (reactant_atm.get_positions() - final_product.get_positions()) ** 2
        rmsd = np.sqrt(np.mean(np.sum(diff_sq, axis=1)))
        print(f"Post-minimization RMSD: {rmsd:.6f} Å")

        ase.io.write(output.reactant, reactant_atm)
        ase.io.write(output.product, final_product)

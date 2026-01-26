# -*- mode:snakemake; -*-


rule do_minimization:
    input:
        config="resources/config_minim.ini",
        endpoint=f"{config['paths']['endpoints']}/{{system}}/{{endpoint}}_pre_aligned.con",
        model=expand(
            f"{config['paths']['models']}/pet-mad-{{version}}.pt",
            version=config["pet_mad"]["version"],
        ),
    output:
        endpoint=f"{config['paths']['endpoints']}/{{system}}/{{endpoint}}_minimized.con",
    shadow:
        "minimal"
    shell:
        """
        cp -f {input.config} config.ini
        cp -f {input.endpoint} pos.con
        cp -f {input.model} .
        eonclient
        cp -f min.con {output.endpoint}
        """


rule do_neb:
    input:
        # Select config based on the {method} wildcard
        config=lambda wildcards: config["methods"][wildcards.method]["neb"],
        reactant=f"{config['paths']['endpoints']}/{{system}}/reactant.con",
        product=f"{config['paths']['endpoints']}/{{system}}/product.con",
        # Use a lambda to access 'wildcards.system'
        path_images=lambda w: expand(
            f"{config['paths']['idpp']}/{w.system}/path/{{num:02d}}.con",
            num=range(config["common"]["number_of_intermediate_imgs"] + 2),
        ),
        model=expand(
            f"{config['paths']['models']}/pet-mad-{{version}}.pt",
            version=config["pet_mad"]["version"],
        ),
    output:
        results_dat=f"{config['paths']['neb']}/{{system}}/{{method}}/results.dat",
        neb_con=f"{config['paths']['neb']}/{{system}}/{{method}}/neb.con",
        neb_dat=f"{config['paths']['neb']}/{{system}}/{{method}}/neb.dat",
    params:
        opath=f"{config['paths']['neb']}/{{system}}/{{method}}",
    shell:
        """
        rm -rf {params.opath}
        mkdir -p {params.opath}
        cp {input.model} {params.opath}/
        cp {input.config} {params.opath}/config.ini
        cp {input.reactant} {params.opath}/reactant.con
        cp {input.product} {params.opath}/product.con

        cd {params.opath}
        eonclient 2>&1 || true
        """


rule do_dimer:
    input:
        config=lambda w: config["methods"][w.method]["mmf"],
        # Depend on NEB completion
        neb_result=f"{config['paths']['neb']}/{{system}}/{{method}}/results.dat",
        model=expand(
            f"{config['paths']['models']}/pet-mad-{{version}}.pt",
            version=config["pet_mad"]["version"],
        ),
    output:
        results_dat=f"{config['paths']['mmf']}/{{system}}/{{method}}/results.dat",
        saddle_con=f"{config['paths']['mmf']}/{{system}}/{{method}}/saddle.con",
    params:
        opath=f"{config['paths']['mmf']}/{{system}}/{{method}}",
        nebroot=f"{config['paths']['neb']}/{{system}}/{{method}}",
    run:
        from ase.io import read, write
        import numpy as np
        import os
        import glob
        import re
        import subprocess
        import shutil

        target_dir = params.opath
        if os.path.exists(target_dir):
            shutil.rmtree(target_dir)
        os.makedirs(target_dir, exist_ok=True)

        for m in input.model:
            shutil.copy(m, target_dir)
        shutil.copy(input.config, os.path.join(target_dir, "config.ini"))

        # Locate peak files in the NEB output directory
        neb_dir = params.nebroot
        mode_files = glob.glob(os.path.join(neb_dir, "peak*_mode.dat"))
        mode_files.sort()

        if not mode_files:
            raise Exception(
                f"No peak*_mode.dat files found in {neb_dir}. Cannot seed Dimer."
            )

            # Use the first peak (peak00)
        mode_file = mode_files[0]
        match = re.search(r"peak(\d+)_mode\.dat", os.path.basename(mode_file))
        if not match:
            raise Exception(f"Failed to parse index from {mode_file}")

        idx = match.group(1)
        orig_pos = os.path.join(neb_dir, f"peak{idx}_pos.con")

        if not os.path.exists(orig_pos):
            raise Exception(f"Missing position file {orig_pos} for mode {mode_file}")

        print(f"[{wildcards.system}] Seeding Dimer from NEB Peak {idx}...")

        # a. Create pos.con
        atoms = read(orig_pos)
        write(os.path.join(target_dir, "pos.con"), atoms)

        # b. Create direction.dat
        direction_vector = np.loadtxt(mode_file)
        np.savetxt(
            os.path.join(target_dir, "direction.dat"), direction_vector, fmt="%12.8f"
        )

        # c. Create displacement.con (finite difference)
        delta = 1e-2
        displaced_atoms = atoms.copy()
        displaced_atoms.positions += delta * direction_vector
        write(os.path.join(target_dir, "displacement.con"), displaced_atoms)

        # 3. Run Dimer
        subprocess.run(["eonclient"], cwd=target_dir, check=True)

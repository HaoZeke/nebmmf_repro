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
        config=lambda wildcards: config["methods"][wildcards.method],
        idpp_path=f"{config['paths']['idpp']}/{{system}}/idppPath.dat",
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
        cp {input.idpp_path} {params.opath}/
        cp {input.reactant} {params.opath}/reactant.con
        cp {input.product} {params.opath}/product.con

        cd {params.opath}
        eonclient 2>&1 || true
        """

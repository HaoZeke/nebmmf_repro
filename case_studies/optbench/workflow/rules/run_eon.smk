# -*- mode:snakemake; -*-


rule do_neb:
    """
    Executes the EON client to find the minimum energy path (MEP).
    This rule pulls source geometries directly from the config to avoid staging.
    """
    input:
        # Fetch method-specific config and system-specific geometries directly
        config_ini=lambda wildcards: config["methods"][wildcards.method],
        reactant_src=lambda wildcards: config["systems"][wildcards.system]["reactant"],
        product_src=lambda wildcards: config["systems"][wildcards.system]["product"],
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
        
        # Copy original files directly into the execution path
        cp {input.config_ini} {params.opath}/config.ini
        cp {input.reactant_src} {params.opath}/reactant.con
        cp {input.product_src} {params.opath}/product.con

        cd {params.opath}
        eonclient 2>&1 || true
        """

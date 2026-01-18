# -*- mode:snakemake; -*-


# 1. Aggregation Rule
rule all_neb_plots:
    input:
        plot1d_path=f"{config['paths']['plots']}/1D_path.png",
        plot1d_index=f"{config['paths']['plots']}/1D_index.png",
    output:
        touch(f"{config['paths']['plots']}/.neb_plots.done"),


# 2. Recipe for 1D Path (Reaction Coordinate: Path)
rule plot_neb_1d_path:
    input:
        con=f"{config['paths']['neb']}/{{system}}/{{method}}/neb.con",
    output:
        plot=f"{config['paths']['plots']}/{{system}}/{{method}}/1D_path.png",
    params:
        ipath=f"{config['paths']['neb']}/{{system}}/{{method}}",
        _cache=f"{config['paths']['cache']}/{{system}}/{{method}}",
        plot_structures=config["plotting"]["plot_structures"],
        facecolor=config["plotting"]["facecolor"],
        figsize=config["plotting"]["figure"]["figsize"],
        dpi=config["plotting"]["figure"]["dpi"],
        fontsize_base=config["plotting"]["fonts"]["base"],
        title=lambda w: f"{config['plotting']['title']}{w.method}",
        zoom_ratio=config["plotting"]["figure"]["zoom_ratio"],
        # Geometry drawings
        dp_x=config["plotting"]["draw_product"]["x"],
        dp_y=config["plotting"]["draw_product"]["y"],
        dp_rad=config["plotting"]["draw_product"]["rad"],
        dr_x=config["plotting"]["draw_reactant"]["x"],
        dr_y=config["plotting"]["draw_reactant"]["y"],
        dr_rad=config["plotting"]["draw_reactant"]["rad"],
        ds_x=config["plotting"]["draw_saddle"]["x"],
        ds_y=config["plotting"]["draw_saddle"]["y"],
        ds_rad=config["plotting"]["draw_saddle"]["rad"],
        # Misc
        _aserot=config["plotting"]["aserot"],
    shell:
        """
        mkdir -p {params._cache} &&
        python -m rgpycrumbs.cli eon plt_neb \
            --con-file {input.con} \
            --output-file {output.plot} \
            --plot-type "profile" \
            --rc-mode "path" \
            --plot-structures "{params.plot_structures}" \
            --facecolor "{params.facecolor}" \
            --input-dat-pattern "{params.ipath}/neb_*.dat" \
            --figsize {params.figsize} \
            --dpi {params.dpi} \
            --zoom-ratio {params.zoom_ratio} \
            --fontsize-base {params.fontsize_base} \
            --draw-product {params.dp_x} {params.dp_y} {params.dp_rad} \
            --draw-reactant {params.dr_x} {params.dr_y} {params.dr_rad} \
            --draw-saddle {params.ds_x} {params.ds_y} {params.ds_rad} \
            --title "{params.title}" \
            --ase-rotation {params._aserot}
        """


# 3. Recipe for 1D Index (Reaction Coordinate: Image Index)
rule plot_neb_1d_index:
    input:
        con=f"{config['paths']['neb']}/{{system}}/{{method}}/neb.con",
    output:
        plot=f"{config['paths']['plots']}/{{system}}/{{method}}/1D_index.png",
    params:
        ipath=f"{config['paths']['neb']}/{{system}}/{{method}}",
        _cache=f"{config['paths']['cache']}/{{system}}/{{method}}",
        plot_structures=config["plotting"]["plot_structures"],
        facecolor=config["plotting"]["facecolor"],
        figsize=config["plotting"]["figure"]["figsize"],
        dpi=config["plotting"]["figure"]["dpi"],
        fontsize_base=config["plotting"]["fonts"]["base"],
        title=lambda w: f"{config['plotting']['title']}{w.method}",
        zoom_ratio=config["plotting"]["figure"]["zoom_ratio"],
        dp_x=config["plotting"]["draw_product"]["x"],
        dp_y=config["plotting"]["draw_product"]["y"],
        dp_rad=config["plotting"]["draw_product"]["rad"],
        dr_x=config["plotting"]["draw_reactant"]["x"],
        dr_y=config["plotting"]["draw_reactant"]["y"],
        dr_rad=config["plotting"]["draw_reactant"]["rad"],
        ds_x=config["plotting"]["draw_saddle"]["x"],
        ds_y=config["plotting"]["draw_saddle"]["y"],
        ds_rad=config["plotting"]["draw_saddle"]["rad"],
        _aserot=config["plotting"]["aserot"],
    shell:
        """
        mkdir -p {params._cache} &&
        python -m rgpycrumbs.cli eon plt_neb \
            --con-file {input.con} \
            --output-file {output.plot} \
            --plot-type "profile" \
            --rc-mode "index" \
            --plot-structures "{params.plot_structures}" \
            --facecolor "{params.facecolor}" \
            --input-dat-pattern "{params.ipath}/neb_*.dat" \
            --figsize {params.figsize} \
            --dpi {params.dpi} \
            --zoom-ratio {params.zoom_ratio} \
            --fontsize-base {params.fontsize_base} \
            --draw-product {params.dp_x} {params.dp_y} {params.dp_rad} \
            --draw-reactant {params.dr_x} {params.dr_y} {params.dr_rad} \
            --draw-saddle {params.ds_x} {params.ds_y} {params.ds_rad} \
            --title "{params.title}" \
            --ase-rotation {params._aserot}
        """

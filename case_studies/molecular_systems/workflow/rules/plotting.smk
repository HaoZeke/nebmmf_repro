# -*- mode:snakemake; -*-


# 1. Aggregation Rule
rule all_neb_plots:
    input:
        plot1d_path=f"{config['paths']['plots']}/1D_path.png",
        plot1d_index=f"{config['paths']['plots']}/1D_index.png",
        plot1d_rmsd=f"{config['paths']['plots']}/1D_rmsd.png",
        plot2d_rmsd=f"{config['paths']['plots']}/2D_rmsd.png",
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


# 4. Recipe for 1D RMSD (Reaction Coordinate: RMSD)
rule plot_neb_1d_rmsd:
    input:
        con=f"{config['paths']['neb']}/{{system}}/{{method}}/neb.con",
    output:
        plot=f"{config['paths']['plots']}/{{system}}/{{method}}/1D_rmsd.png",
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
            --rc-mode "rmsd" \
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
            --cache-file {params._cache}/1dcache.parquet \
            --ase-rotation {params._aserot}
        """


# 5. Recipe for 2D Landscape (Reaction Coordinate: RMSD)
# Can't fail the output rules
# To add the reference structures use
# --additional-con {params._addl} "baker" \
rule plot_neb_2d_rmsd:
    input:
        con=f"{config['paths']['neb']}/{{system}}/{{method}}/neb.con",
    output:
        plot=f"{config['paths']['plots']}/{{system}}/{{method}}/2D_rmsd.png",
    params:
        ipath=f"{config['paths']['neb']}/{{system}}/{{method}}",
        _cache=f"{config['paths']['cache']}/{{system}}/{{method}}",
        plot_structures=config["plotting"]["plot_structures"],
        facecolor=config["plotting"]["facecolor"],
        figsize=config["plotting"]["figure"]["figsize"],
        dpi=config["plotting"]["figure"]["dpi"],
        fontsize_base=config["plotting"]["fonts"]["base"],
        title=lambda w: f"RMSD(R,P) projection ({w.method})",
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
        _addl=lambda w: config["systems"][w.system]["saddle"],
        _additional_flags=lambda w: " ".join([
            f'--additional-con {config["paths"]["neb"]}/{w.system}/{m}/sp.con "{m.upper()}"'
            for m in ["cineb", "mmf"] if m != w.method
        ])
    shell:
        """
        mkdir -p {params._cache} &&
        python -m rgpycrumbs.cli eon plt_neb \
            --con-file {input.con} \
            --output-file {output.plot} \
            --plot-type "landscape" \
            --rc-mode "rmsd" \
            --show-pts \
            --landscape-path "all" \
            --plot-structures "{params.plot_structures}" \
            --facecolor "{params.facecolor}" \
            --input-dat-pattern "{params.ipath}/neb_*.dat" \
            --input-path-pattern "{params.ipath}/neb_path*.con" \
            --figsize {params.figsize} \
            --dpi {params.dpi} \
            --zoom-ratio {params.zoom_ratio} \
            --fontsize-base {params.fontsize_base} \
            --draw-product {params.dp_x} {params.dp_y} {params.dp_rad} \
            --draw-reactant {params.dr_x} -{params.dr_y} {params.dr_rad} \
            --draw-saddle {params.ds_x} {params.ds_y} {params.ds_rad} \
            --title "{params.title}" \
            --cache-file {params._cache}/2dcache.parquet \
            --ase-rotation {params._aserot} \
            --ira-kmax 14 \
            {params._additional_flags} \
            --show-legend || touch {output.plot}
        """

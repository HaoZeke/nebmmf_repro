# -*- mode:snakemake; -*-
"""
Post-NEB analysis rules: parse results, brms Bayesian analysis,
summary figures, case studies, static handover, archive.
"""

from pathlib import Path

NEB_DIR = config["paths"]["neb"]
BAKER_SYSTEMS = list(config.get("systems", {}).keys()) if isinstance(
    config.get("systems"), dict
) else config.get("systems", [])
METHODS = config.get("methods", [])
PAPER_IMGS = "/home/goswami/Git/Github/TeX/roneb_tex/v2/FrontiersChemistry/imgs"


# ---------------------------------------------------------------------------
# 1. Parse NEB results into benchmark CSV (with IRA RMSD)
# ---------------------------------------------------------------------------
rule parse_results:
    input:
        expand(f"{NEB_DIR}/{{sys}}/{{method}}/results.dat",
               sys=BAKER_SYSTEMS, method=METHODS),
    output:
        csv="data/baker_bench.csv",
    shell:
        "python ../scripts/parse_results.py {NEB_DIR} -o {output.csv}"


# ---------------------------------------------------------------------------
# 2. Optuna parameter sensitivity (long-running, NOT in rule all)
# ---------------------------------------------------------------------------
rule optuna_study:
    """Run separately: pixi run -e eongpu snakemake optuna_study"""
    output:
        db="../results/optuna_study.db",
        params_yaml="../results/optuna_best_params.yaml",
    shell:
        "python ../scripts/optuna_param_study.py --n-trials 200"

rule optuna_2param:
    """Run separately: pixi run -e eongpu snakemake optuna_2param"""
    output:
        db="../results/optuna_2param.db",
    shell:
        "python ../scripts/optuna_2param.py --n-trials 200"


# ---------------------------------------------------------------------------
# 3. brms Bayesian reanalysis (R, brms pixi env)
# ---------------------------------------------------------------------------
rule brms_analysis:
    input:
        csv="data/baker_bench.csv",
    output:
        effects="data/models/effect_summary.txt",
        model_output="data/models/brms_model_output.txt",
        brms_pes="imgs/gen/R/brms_pes.png",
        pp_density="imgs/gen/R/suppl/pp_density.png",
        pp_group="imgs/gen/R/suppl/pp_group.png",
        shape_post="imgs/gen/R/suppl/brms_shape_posterior.png",
        ppc_loo="imgs/gen/R/suppl/ppc_loo.png",
        nolog_dist="imgs/gen/R/suppl/nolog_dist.png",
        violin="imgs/gen/R/suppl/cactus_violin.png",
    shell:
        "pixi run -e brms Rscript --no-init-file ../scripts/brms_reanalysis.R"


# ---------------------------------------------------------------------------
# 4. Dumbbell comparison plot (R, rviz env)
# ---------------------------------------------------------------------------
rule dumbbell_plot:
    input:
        csv="data/baker_bench.csv",
    output:
        plot="imgs/gen/R/dumbbell_plot.png",
    shell:
        "pixi run -e rviz Rscript --no-init-file ../scripts/regen_dumbbell.R"


# ---------------------------------------------------------------------------
# 5. Dataset characterization (R, brms env)
# ---------------------------------------------------------------------------
rule dataset_characterization:
    input:
        csv="data/baker_bench.csv",
    output:
        plot="imgs/gen/R/suppl/dataset_characterization.png",
    shell:
        "pixi run -e brms Rscript --no-init-file ../scripts/gen_dataset_char.R"


# ---------------------------------------------------------------------------
# 6. Supplementary baker_2d (48 2D landscapes + 48 1D profiles)
# ---------------------------------------------------------------------------
rule supplementary_figures:
    input:
        expand(f"{NEB_DIR}/{{sys}}/{{method}}/neb.con",
               sys=BAKER_SYSTEMS, method=METHODS),
        csv="data/baker_bench.csv",
    output:
        touch("results/.baker_2d_done"),
    shell:
        "python ../scripts/gen_baker_2d.py"


# ---------------------------------------------------------------------------
# 7. Static handover comparison (Claisen System 17)
# ---------------------------------------------------------------------------
rule static_handover:
    input:
        reactant=f"{config['paths']['endpoints']}/17_claisen/reactant.con",
        product=f"{config['paths']['endpoints']}/17_claisen/product.con",
    output:
        touch("results/static_handover/.done"),
    shell:
        "python ../scripts/run_static_handover.py"


# ---------------------------------------------------------------------------
# 8. Case study plots (HNCCS + Claisen + HCONH3+, main + suppl figures)
# ---------------------------------------------------------------------------
rule case_study_plots:
    input:
        expand(f"{NEB_DIR}/{{sys}}/{{method}}/neb.con",
               sys=["01_hcn", "19_hnccs", "17_claisen", "20_hconh3_cation"], method=METHODS),
    output:
        touch("results/.case_study_done"),
    run:
        import subprocess, os, tempfile
        BASE_CWD = os.getcwd()
        # Each entry: (sys, method, ptype, outfile, title, [(con_path, label), ...])
        plots = [
            # HNCCS (main paper Figure 5)
            ("19_hnccs", "mmf", "landscape", f"{PAPER_IMGS}/hnccs_mmf_2D.png",
             "19: HNCCS (OCI-NEB)",
             [(f"{NEB_DIR}/19_hnccs/cineb/sp.con", "CINEB")]),
            ("19_hnccs", "cineb", "profile", f"{PAPER_IMGS}/hnccs_cineb_1D_path.png",
             "19: HNCCS (CI-NEB)", []),
            ("19_hnccs", "mmf", "profile", f"{PAPER_IMGS}/hnccs_mmf_1D_path.png",
             "19: HNCCS (OCI-NEB)", []),
            # HCN static failure case study (main paper): dimer finds WRONG saddle
            # OCI marker omitted from strip (overlaps SP); both shown on landscape
            ("01_hcn", "cineb", "landscape", f"{PAPER_IMGS}/hcn_static_failure_2D.png",
             "01: HCN",
             [("results/static_handover/01_hcn/dimer/saddle.con", "Static")]),
            # Claisen static comparison (supplementary, both methods converge)
            ("17_claisen", "cineb", "landscape", f"{PAPER_IMGS}/cases/17_claisen_static_neb_2D.png",
             "17: Claisen",
             [(f"{NEB_DIR}/17_claisen/mmf/sp.con", "OCINEB"),
              ("results/static_handover/17_claisen/dimer/saddle.con", "Dimer")]),
            # HCONH3+ cation ablation (supplementary)
            ("20_hconh3_cation", "cineb", "landscape", f"{PAPER_IMGS}/cases/hconh3cat_cineb_2D.png",
             "20: HCONH3+ (CI-NEB)", []),
            ("20_hconh3_cation", "mmf", "profile", f"{PAPER_IMGS}/cases/hconh3cat_mmf_agg_1D_path.png",
             "20: HCONH3+ (OCI-NEB)", []),
            ("20_hconh3_cation", "cineb", "profile", f"{PAPER_IMGS}/cases/hconh3cat_mmf_strict_1D_path.png",
             "20: HCONH3+ (CI-NEB)", []),
        ]
        # Per-plot strip spacing: more points need more room
        strip_overrides = {"hcn_static_failure_2D.png": "3.5"}
        for sys_name, method, ptype, outfile, title, add_cons in plots:
            os.makedirs(os.path.dirname(outfile), exist_ok=True)
            # Resolve all paths to absolute so they work from a tmpdir
            con_file = os.path.join(BASE_CWD, f"{NEB_DIR}/{sys_name}/{method}/neb.con")
            dat_pat = os.path.join(BASE_CWD, f"{NEB_DIR}/{sys_name}/{method}/neb_*.dat")
            out_basename = os.path.basename(outfile)
            spacing = strip_overrides.get(out_basename, "2.5")
            cmd = [
                "python", "-m", "rgpycrumbs.cli", "--dev", "eon", "plt-neb",
                "--con-file", con_file,
                "--output-file", outfile,
                "--plot-type", ptype,
                "--rc-mode", "rmsd" if ptype == "landscape" else "path",
                "--input-dat-pattern", dat_pat,
                "--plot-structures", "crit_points",
                "--dpi", "300", "--rotation", "auto",
                "--strip-renderer", "xyzrender", "--title", title,
                "--strip-dividers", "--strip-spacing", spacing,
            ]
            if ptype == "landscape":
                path_pat = os.path.join(BASE_CWD, f"{NEB_DIR}/{sys_name}/{method}/neb_path*.con")
                cmd.extend([
                    "--surface-type", "grad_imq", "--show-pts",
                    "--landscape-path", "all",
                    "--input-path-pattern", path_pat,
                    "--ira-kmax", "14", "--show-legend",
                    "--figsize", "5.37", "5.37", "--zoom-ratio", "0.2",
                ])
                for con_path, con_label in add_cons:
                    abs_con = os.path.join(BASE_CWD, con_path) if not os.path.isabs(con_path) else con_path
                    if os.path.exists(abs_con):
                        cmd.extend(["--additional-con", abs_con, con_label])
            else:
                cmd.extend(["--figsize", "5.37", "3.37", "--zoom-ratio", "0.1", "--no-legend"])
            # Run in per-plot tmpdir to isolate .neb_landscape.parquet cache
            with tempfile.TemporaryDirectory(prefix=f"case_{sys_name}_{method}_") as tmpdir:
                subprocess.run(cmd, cwd=tmpdir, check=True, timeout=300)


# ---------------------------------------------------------------------------
# 9. Copy generated R plots to paper directory
# ---------------------------------------------------------------------------
rule sync_paper_figures:
    input:
        brms_pes="imgs/gen/R/brms_pes.png",
        dumbbell="imgs/gen/R/dumbbell_plot.png",
        dataset="imgs/gen/R/suppl/dataset_characterization.png",
        violin="imgs/gen/R/suppl/cactus_violin.png",
        pp_density="imgs/gen/R/suppl/pp_density.png",
        pp_group="imgs/gen/R/suppl/pp_group.png",
        shape_post="imgs/gen/R/suppl/brms_shape_posterior.png",
        ppc_loo="imgs/gen/R/suppl/ppc_loo.png",
        nolog_dist="imgs/gen/R/suppl/nolog_dist.png",
    output:
        touch("results/.paper_sync_done"),
    shell:
        """
        mkdir -p {PAPER_IMGS}/suppl
        cp {input.brms_pes} {PAPER_IMGS}/brms_pes.png
        cp {input.dumbbell} {PAPER_IMGS}/dumbbell_plot.png
        cp {input.dataset} {PAPER_IMGS}/dataset_characterization.png
        cp {input.violin} {PAPER_IMGS}/suppl/cactus_violin.png
        cp {input.pp_density} {PAPER_IMGS}/suppl/pp_density.png
        cp {input.pp_group} {PAPER_IMGS}/suppl/pp_group.png
        cp {input.shape_post} {PAPER_IMGS}/suppl/brms_shape_posterior.png
        cp {input.ppc_loo} {PAPER_IMGS}/suppl/ppc_loo.png
        cp {input.nolog_dist} {PAPER_IMGS}/suppl/nolog_dist.png
        """


# ---------------------------------------------------------------------------
# 10. Materials Cloud archive
# ---------------------------------------------------------------------------
ARCHIVE = "results/nebmmf_archive.tar.xz"

rule archive:
    input:
        csv="data/baker_bench.csv",
        brms="imgs/gen/R/brms_pes.png",
        dumbbell="imgs/gen/R/dumbbell_plot.png",
        baker2d="results/.baker_2d_done",
        cases="results/.case_study_done",
        static="results/static_handover/.done",
        dataset="imgs/gen/R/suppl/dataset_characterization.png",
        paper_sync="results/.paper_sync_done",
    output:
        ARCHIVE,
    shell:
        """
        tar -cJf {output} \
            --exclude='.snakemake' --exclude='__pycache__' \
            --exclude='*.egg-info' --exclude='nebmmf_archive.tar.xz' \
            data/ results/ config/ workflow/ imgs/gen/ \
            ../scripts/ ../results/ \
            ../pixi.toml ../pixi.lock
        """

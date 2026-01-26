# Preliminary

Clone the GitHub repository:

``` bash
git clone git@github.com:HaoZeke/nebmmf_repro.git
```

# About

Contains the reproduction details for the publication on the performance and
logs for the RONEB algorithm / OCI-NEB for seeking accelerated transition state.

## Reference

If you use this repository or its parts please cite the corresponding publication or data source.


### Preprint

> 1.  Goswami, M. Gunde, and H. Jónsson, “Enhanced climbing image nudged elastic band method with hessian eigenmode alignment,” Jan. 22, 2026, arXiv: arXiv:2601.12630. doi: 10.48550/arXiv.2601.12630.


## Replication data

Users must inflate the archives from the MaterialsCloud record into the
repository structure. Place the following archives into the root or the
designated subdirectories before execution:

  Archive                          Destination                        Description
  -------------------------------- ---------------------------------- ----------------------------------------------------------------
  `data.tar.xz`                    `data/`                            Stan models, Parquet data, and CSV benchmarks.
  `case_baker_results.tar.xz`      `case_studies/molecular_systems`   Parameter ablation results for the Baker benchmark set.
  `optbench_results.tar.xz`        `case_studies/optbench`            Platinum heptamer benchmark for surface systems from OptBench.
  `case_static_results.tar.xz`     `case_studies/static_switchover`   Baker test set, one-shot NEB followed by Dimer.
  `eonRuns_results.tar.xz`         `eonRuns/`                         Full NEB/RONEB trajectory logs and plots for the Baker.

Inflate the results archive using `tar`, e.g.

``` bash
tar -xvf eonRuns_results.tar.xz -C eonRuns/
```

### Reusing models

The repository offers both F.A.I.R.-compliant artifacts and native `R` objects for Bayesian analysis via `brms`.

1.  From R

    The `data/models/` directory stores serialized `.rds` objects,
    readable with `readRDS`. These require the `brms` library for
    further analysis:

    ``` {.r org-language="R" eval="never"}
    library('brms')
    model <- readRDS("data/models/brms_efficiency_scaling_v5.rds")
    ```

    More helper functions for generating and using these models and
    predictions are in the Github repository.

    The repository employs `pixi` for dependency management and
    organizes content by function:

      Directory         Description
      ----------------- --------------------------------------------------------------------------------
      `case_studies/`   Specific molecular systems and `optbench` configurations.
      `data/`           Raw benchmarks (`.csv`), Stan models, and Parquet data.
      `docs/`           Documentation, including `.org` files for visualization and model definitions.
      `eonRuns/`        Workflow configurations and results for the eOn runs.
      `scripts/`        Python utilities for generating YAML configs and parsing simulation outputs.

    The results directories (`eonRuns/results`) implement the following hierarchy:

    -   ****00~models~****: Contains the PyTorch (`.pt`) Machine
        Learning Potentials used for the PES evaluations.
    -   ****01~endpoints~****: Atomic coordinates for the reactant and
        product states for each system (e.g., `01_hcn`,
        `09_parentdielsalder`).
    -   ****02~idpp~****: Initial guesses for the reaction paths.
    -   ****03~neb~****: Optimized transition state paths using the
        RONEB and standard NEB algorithms.
    -   ****04~plots~****: Resultant visualizations and energy profiles
        for each chemical system.

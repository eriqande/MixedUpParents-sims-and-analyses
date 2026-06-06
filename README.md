
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MixedUpParents-sims-and-analyses

The full workflow runs the large simulation studies presented in the
paper “MixedUpParents: A New and Improved Method for Genetic Parentage
Assignment in Admixed Populations.” For a smaller reproducibility check,
the Snakefile includes a `README_example` rule. It requests four
ROC-curve outputs from replicate 2 of one `cyclone_nonWF` scenario:

- `naive_logl_both_diag_and_var_rocs.rds`
- `mup_rocs.rds`
- `hot_both_diag_and_var_rocs.rds`
- `relate_admix_both_diag_and_var_rocs.rds`

Run a dry-run of that example with:

``` sh
snakemake -np --use-conda  --cores 8 \
  README_example \
  --config mup_conda="$HOME/miniforge3/envs/mixed-up-parents"
```

Remove `-n` to actually run the example, which will take about 30
minutes on a typical laptop.

The `mup_conda` config value should point to the conda environment used
by `MixedUpSlimSims`/`MixedUpParents`; the workflow exports it as
`MUP_CONDA` for `reticulate`.

The RelateAdmix part of the example also expects `plink`, `admixture`,
and `relateAdmix` to be available on `PATH`. Those dependencies are not
handled via conda.

The MUP log-likelihood rule uses Snakemake’s `--use-conda` machinery to
create the small bedtools environment in `envs/bedtools.yaml`.

## Full dry run on main branch

A full dry-run can be done on the main branch with:

``` sh
snakemake -np --config mup_conda=$HOME/miniforge3/envs/mixed-up-parents
```

Again, the mup_conda config must point to the `mixed-up-parents` conda
env.

Such a run should produce:

    Job stats:
    job             count
    ------------  -------
    all                 1
    compute_rocs     3200
    gather_rocs         1
    hot_scores       1280
    mup_logls         640
    naive_logls      1280
    run_sequoia       300
    slim_sim           40
    tweak2mup         640
    total            7382

The full run uses `config/config.yaml`, which, upon inspection, you can
see references `config/sim_spec.csv`. This is a CSV-based table that
dictates how many replicates of each particular combination of
population size, migration rate, sampled proportion, genotyping error
rate, and missing genotype rate is to be simulated. The first ten lines
of that file are:

``` sh
sim_scenario,num_reps,pop_size_1,pop_size_2,mig_rate_1,mig_rate_2,prop_sampled,err_var,err_diag,miss_var,miss_diag
cyclone_WF,5,1200,1200,0,0,1,0.01,0.004,0,0
cyclone_WF,5,1200,1200,0,0,1,0.01,0.004,0.15,0.15
cyclone_WF,5,1200,1200,0,0,1,0.01,0.004,0.25,0.25
cyclone_WF,5,1200,1200,0,0,1,0.01,0.004,0.5,0.5
cyclone_WF,5,1200,1200,0,0,0.75,0.01,0.004,0,0
cyclone_WF,5,1200,1200,0,0,0.75,0.01,0.004,0.15,0.15
cyclone_WF,5,1200,1200,0,0,0.75,0.01,0.004,0.25,0.25
cyclone_WF,5,1200,1200,0,0,0.75,0.01,0.004,0.5,0.5
cyclone_WF,5,1200,1200,0,0,0.5,0.01,0.004,0,0
```

See the section `Making a sim-spec.csv file` for information on how to
assemble such a file.

## Separate branches

There are two additional branches of this
repository—`rerun-sequoia-smalls` and `relate-admix-on-alpine`—that hold
configs for doing additional simulations or running additional inference
routines, respectively.

`rerun-sequoia-smalls` is maintained for running the small-sample cases,
and `relate-admix-on-alpine` was used to run RelateAdmix on the Alpine
cluster.

## Plotting Tools

The ROC figures appearing in the paper were made from the gathered
parentage ROC curves produced by rule `gather_rocs` in the notebooks
`004-ROCs-for-big-runs.Rmd` for the full-scale runs and
`004.1--ROCs-for-the-small-120-runs.Rmd` for the smaller set of runs
with a higher chance of Sequoia finishing in less than 24 hours.

The run time boxplots in the paper were created within the notebook
`007-harvest-and-present-benchmarks.Rmd` using Snakemake’s benchmark
outputs for each rule invocation.

## Notes re: Running on the SEDNA Cluster

Something weird goes on with SEDNA whereby the slim_sim rule writes all
its output no problem, and then fails when returning control back to
Snakemake. It might have something to do with the reticulate package and
its conda environment.

At any rate, you can run it like this (which is assuming that we are
running snakemake on a node with 10 cores, 9 of which we will use for
local jobs):

``` sh
snakemake -np --use-envmodules --keep-incomplete --until slim_sim --profile hpcc-profiles/slurm/sedna
```

and once that is done, the outputs will still be there and you can
continue by simply removing the `.snakemake/incomplete` directory and
then doing:

``` sh
snakemake -np --use-envmodules --use-conda --local-cores 19 --profile hpcc-profiles/slurm/sedna
```

Total PITA. Gonna have to figure out why it is failing, but we can use
the above until that time.

We might be able to use –cleanup-metadata

Note, the –use-conda on the second line is to have bedtools for
mup_logls.

You also need to pass the `mup_conda` config value, as in the README
example above.

## Making a sim-spec.csv file

This is just some example code on how to make these

``` r
library(tidyverse)
sim_specs <- expand_grid(
  sim_scenario = c("cyclone_WF", "cyclone_nonWF"),
  num_reps = 5,
  pop_size_1 = 1200,
  pop_size_2 = 1200,
  nesting(
    mig_rate_1 = c(0.0, 0.02, 0.05, 0.1),
    mig_rate_2 = c(0.0, 0.02, 0.05, 0.1),
  ),
  prop_sampled = c(1.0, 0.75, 0.5, 0.25),
  err_var = 0.01,
  err_diag = 0.004,
  nesting(
    miss_var = c(0, 0.15, 0.25, 0.5),
    miss_diag = c(0, 0.15, 0.25, 0.5),
  )
)
```

## Making the small (120 indivs) simulations

Only do certain Sequoia runs. Based on previous runs the only ones that
will finish in a reasonable amount of time appear to be:

- exclude_same_cohort
- var only
- missing locus proportion of 0, 0.15

I hope that will finish. This is all specified in a separate config
called `config/config_120.yaml`.

``` r
library(tidyverse)
sim_specs <- expand_grid(
  sim_scenario = c("small_cyclone_WF", "small_cyclone_nonWF"),
  num_reps = 5,
  pop_size_1 = 120,
  pop_size_2 = 120,
  nesting(
    mig_rate_1 = c(0.0, 0.02, 0.05, 0.1),
    mig_rate_2 = c(0.0, 0.02, 0.05, 0.1),
  ),
  prop_sampled = c(0.5),
  err_var = 0.01,
  err_diag = 0.004,
  nesting(
    miss_var = c(0.25, 0.5),
    miss_diag = c(0.25, 0.5),
  )
)
```

Even at that, with full sampling of individuals we might have some
really long ones. But at least that will give us only 2 \* 5 \* 4 \* 4
\* 2 = 320 runs of sequoia, and they should be shorter

I will modify the aggregation functions to have an option to do only
those.

I did that, and I was able to run it with:

``` sh
snakemake -np --use-envmodules --use-conda --local-cores 3  --configfile config/config_120.yaml  --profile hpcc-profiles/slurm/sedna
```

5 of the sequoia runs had not finished in 48 hours. So I just removed
them from the targets. Those files are in

    config/black-lists/sequoia-short-run-timeouts.txt

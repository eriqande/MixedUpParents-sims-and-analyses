
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MixedUpParents-sims-and-analyses

## Notes on SEDNA

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
continue with:

``` sh
snakemake -np --use-envmodules --use-conda --ignore-incomplete --local-cores 9 --profile hpcc-profiles/slurm/sedna
```

Total PITA. Gonna have to figure out why it is failing, but we can use
the above until that time.

We might be able to use –cleanup-metadata

Note, the –use-conda on the second line is to have bedtools for
mup_logls.

You gotta define a MUP_CONDA variable, too. see the readme for
MixedUpSlimSims for that.

## Notes on scenarios:

**Langford simulations**

Migration rates: symmetrical at 0, 2, 5, and 10% Missing indiv rate: 0,
25, 50, 75% Missing locus rate: 0, 0.15, 0.25, 0.5

5 reps of each.

Don’t bother with Sequoia for the same cohort included runs, none of
them finished!

Calculate number of markers shared for each pair.

**Sequoia Runs with smaller sample sizes**

**Other cases**

Fst = 0, hack things up so that we can tease out how much performace is
due to SOS model vs better-estimated allele freqs.

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

Only do certain Sequoia runs. Based on previous runs the only ones that
will finish in a reasonable amount of time appear to be:

- exclude_same_cohort
- var only
- missing locus proportion of 0, 0.15

Even at that, with full sampling of individuals we might have some
really long ones. But at least that will give us only 2 \* 5 \* 4 \* 4
\* 2 = 320 runs of sequoia, and they should be shorter

I will modify the aggregation functions to have an option to do only
those.

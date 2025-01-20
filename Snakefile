# These rules run pretty quickly, so we make them local rules
localrules: hot_scores,compute_rocs,gather_rocs

# this is a snakefile to orchestrate and run the simulations for the MUP paper.
# We have designed simulations with data that look like trout populations.
# Each simulated data set is analyzed by MixedUpParents, HipHop, and Sequoia.

import pandas as pd

# The runs are configured by using tabular configuration, with that file
# specified in the config file.
configfile: "config/config.yaml"

sim_spec = pd.read_csv(config["sim_spec"])

# here is a function to expand file paths.  I am developing it as we go and
# it will change as I get more rules done
def expand_paths():
  x = sim_spec
  return([ expand("results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/mup_rocs.rds",
    slim = x.loc[i, "sim_scenario"],
    ps1 = x.loc[i, "pop_size_1"],
    ps2 = x.loc[i, "pop_size_2"],
    mr1 = x.loc[i, "mig_rate_1"],
    mr2 = x.loc[i, "mig_rate_2"],
    rep = range(x.loc[i, "num_reps"]),
    ppn = x.loc[i, "prop_sampled"],
    verr = x.loc[i, "err_var"],
    derr = x.loc[i, "err_diag"],
    vmiss = x.loc[i, "miss_var"],
    dmiss = x.loc[i, "miss_diag"]) for i in range(x.shape[0])
  ])


# this version lets me request different files within these directories
# with the "what" parameter
def expand_paths_general(spec = sim_spec, what = "mup_rocs.rds"):
  x = spec
  return([ expand("results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/{w}",
    slim = x.loc[i, "sim_scenario"],
    ps1 = x.loc[i, "pop_size_1"],
    ps2 = x.loc[i, "pop_size_2"],
    mr1 = x.loc[i, "mig_rate_1"],
    mr2 = x.loc[i, "mig_rate_2"],
    rep = range(1, x.loc[i, "num_reps"] + 1),
    ppn = x.loc[i, "prop_sampled"],
    verr = x.loc[i, "err_var"],
    derr = x.loc[i, "err_diag"],
    vmiss = x.loc[i, "miss_var"],
    dmiss = x.loc[i, "miss_diag"],
    w = what) for i in range(x.shape[0])
  ])




# here we get a list of the sequoia files we want
SEQUOIA_ALL_FILES = [
  expand(x, sc = ["exclude_same_cohort", "include_same_cohort"], sm = ["only_var", "both_diag_and_var"]) 
  for x in expand_paths_general(what = "sequoia/{sc}-{sm}-rocresults.rds")]


# here we get things only for the big runs at the moment.  Things we hope will finish
filtered_spec = sim_spec[sim_spec["miss_var"] <= 0.15].reset_index(drop=True)
SEQUOIA_FIRST_RUN = [
  expand(x, sc = ["exclude_same_cohort"], sm = ["only_var"]) 
  for x in expand_paths_general(spec = filtered_spec, what = "sequoia/{sc}-{sm}-rocresults.rds")]



rule all:
  input: "results/summarized/all-rocs.rds"

# simulate a SLiM data set.  SLiM and R seeds are set according to output path name
# currently assumes we are pulling samples from the last three generations (9-11)
rule slim_sim:
  input:
    slimfile=lambda wc: config["sim_scenarios"][wc.slim]["slim_base"],
    varfreqs=lambda wc: config["sim_scenarios"][wc.slim]["var_freqs"],
    diagmarkers=lambda wc: config["sim_scenarios"][wc.slim]["diag_markers"],
  params:
    pop_size1="{ps1}",
    pop_size2="{ps2}",
    mig_rate1="{mr1}",
    mig_rate2="{mr2}",
  #conda:
  #  "/Users/eriq/mambaforge-arm64/envs/mup"
  threads: 2 # so I don't run out of memory on my laptop
  output:
    outrds="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/slim-output.rds"
  log:
    log="results/logs/slim_sim/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}.log"
  benchmark:
    "results/benchmarks/slim_sim/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}.bmk"
  conda:
    "/home/eanderson/mambaforge/envs/mixed-up-parents"
  envmodules:
    "R/4.0.3"
  script:
    "scripts/slim_sim.R"




# a rule to apply genotyping error and subsampling (the tweaking) and then
# spit the genotype data out in a MUP-ready format
rule tweak2mup:
  input:
    inrds="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/slim-output.rds"
  params:
    ppn_sampled="{ppn}",
    var_err="{verr}",
    diag_err="{derr}",
    var_miss="{vmiss}",
    diag_miss="{dmiss}"
  output:
    outrds="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/tweaked2mup.rds"
  log:
    log="results/logs/tweak2mup/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}.log"
  benchmark:
    "results/benchmarks/slim_sim/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}.bmk"
  envmodules:
    "R/4.0.3"
  script:
    "scripts/tweak.R"


# a rule to run Sequoia
# note the wildcards:
#  - seq_cohort: cohort excluded or not
#  - seq_mark:  whether including the diagnostic markers or just the variable ones
rule run_sequoia:
  input:
    inrds="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/tweaked2mup.rds",
  output:
    seq_sim_results="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/sequoia/{seq_cohort}-{seq_mark}-simresults.rds",
    seq_roc_results="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/sequoia/{seq_cohort}-{seq_mark}-rocresults.rds",
  params:
    which_markers="{seq_mark}",
    cohort_situation="{seq_cohort}"
  envmodules:
    "R/4.0.3"
  log:
    log="results/logs/run_sequoia/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/{seq_cohort}-{seq_mark}.log"
  benchmark:
    "results/benchmarks/run_sequoia/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/{seq_cohort}-{seq_mark}.bmk"
  script:
    "scripts/sequoia2.R"  


# get the log-likelihoods for all pairs between the last generation and the last 3-generations.
# I am going to get all of those likelihoods, and then I am going to do the different
# comparisons pertaining to whether you can filter out same-age fish or not by post-processing
# all the pairwise comparisons.  That way we don't have to infer all the ancestry tracts
# multiple times, etc.
rule mup_logls:
  input:
    inrds="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/tweaked2mup.rds"
  output:
    outrds="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/mup_all_pairs.rds",
    outnumLoc="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/mup_num_non_missing_loci.rds"
  threads: 8
  log:
    log="results/logs/mup_logls/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}.log"
  benchmark:
    "results/benchmarks/mup_logls/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}.bmk"
  envmodules:
    "R/4.0.3"
  conda:
    "envs/bedtools.yaml"
  script:
    "scripts/mup_logls.R"


# this is like the mup_logls rule, but it counts HOTs.  The wildcard {hot_markers} can
# be either "both_diag_and_var" or "only_var", and that gets sent to the script as a param.
rule hot_scores:
  input:
    inrds="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/tweaked2mup.rds"
  output:
    outrds="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/hot_{hot_markers}_all_pairs.rds"
  threads: 1
  params:
    marker_set = lambda wc: "both" if "both_diag_and_var" in wc.hot_markers else "var" if "only_var" in wc.hot_markers else None
  log:
    log="results/logs/hot_scores/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/{hot_markers}.log"
  benchmark:
    "results/benchmarks/hot_scores/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/{hot_markers}.bmk"
  envmodules:
    "R/4.0.3"
  script:
    "scripts/hiphop_scores.R"

rule compute_rocs:
  input:
    inrds="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/{what}_all_pairs.rds"
  output:
    outrds="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/{what}_rocs.rds"
  log:
    log="results/logs/compute_rocs/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/{what}.log"
  benchmark:
    "results/benchmarks/compute_rocs/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}/ppn-{ppn}-verr-{verr}-derr-{derr}-vmiss-{vmiss}-dmiss-{dmiss}/{what}.bmk"
  envmodules:
    "R/4.0.3"
  script:
    "scripts/compute_rocs.R"






rule gather_rocs:
  input:
    inList=[
      expand_paths_general(what = "mup_rocs.rds"),
     #expand_paths_general(what = "hot_both_diag_and_var_rocs.rds"),
     # expand_paths_general(what = "hot_only_var_rocs.rds"),
      SEQUOIA_FIRST_RUN 
    ]
  output: 
    outrds="results/summarized/all-rocs.rds"
  log:
    log="results/logs/gather_rocs/log.log"
  benchmark:
    "results/benchmarks/gather_rocs/benchmark.bmk"
  envmodules:
    "R/4.0.3"
  script:
    "scripts/gather_rocs.R"

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
  return([ expand("results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/sim-output/rep-{rep}.rds",
    slim = x.loc[i, "sim_scenario"],
    ps1 = x.loc[i, "pop_size_1"],
    ps2 = x.loc[i, "pop_size_2"],
    mr1 = x.loc[i, "mig_rate_1"],
    mr2 = x.loc[i, "mig_rate_2"],
    rep = range(x.loc[i, "num_reps"])) for i in range(x.shape[0])
  ])



rule all:
  input: expand_paths()

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
  output:
    outrds="results/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/sim-output/rep-{rep}.rds"
  log:
    log="results/logs/slim_sim/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}.log"
  benchmark:
    "results/benchmarks/slim_sim/scenario-{slim}/ps1-{ps1}-ps2-{ps2}-mr1-{mr1}-mr2-{mr2}/rep-{rep}.bmk"
  script:
    "scripts/slim_sim.R"
  
    


# now, a rule to tweak a data set


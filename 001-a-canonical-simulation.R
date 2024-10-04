

# This is just here to do a simple canonical simulation with inference
# to make sure that it is all working correctly, and to get it working
# on the cluster.

library(tidyverse)
library(MixedUpSlimSims)
library(MixedUpParents)

#### simulate data using MixedUpSlimSims ####

simmed <- slim_sim_a_dataset(
  AF = cyclone_creek_var_freqs,
  genome_info = mykiss_chroms,
  founder_pop_nums = c(1000, 1000),
  slim_file = system.file(
    "SLiM-models/2-pop-10-gens-vcf-initialize.slim",
    package = "MixedUpSlimSims"
  ),
  slim_seed = 111,
  years_list = list(`1` = 0:10, `2` = 0:10),
  marker_pop_time_list = list(p1 = 9:11, p2 = 9:11),
  diagnostic_markers = wct_rbt_yct_diagnostic_markers
)

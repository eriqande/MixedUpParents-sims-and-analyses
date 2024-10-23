

# This is just here to do a simple canonical simulation with inference
# to make sure that it is all working correctly, and to get it working
# on the cluster.

library(tidyverse)
library(MixedUpSlimSims)
library(MixedUpParents)

library(reticulate)

# do this explicitly so reticulate does not ask if you want to
# make an environment
use_condaenv(Sys.getenv("MUP_CONDA"))

# simulate data using MixedUpSlimSims
simmed <- slim_sim_a_dataset(
  AF = cyclone_creek_var_freqs,
  genome_info = mykiss_chroms,
  founder_pop_nums = c(1000, 1000),
  slim_file = "development/nonWF-slim/polygamy-neg-binom.slim",
  slim_seed = 111,
  years_list = list(`1` = 0:10, `2` = 0:10),
  marker_pop_time_list = list(p1 = 9:11, p2 = 9:11),
  diagnostic_markers = wct_rbt_yct_diagnostic_markers
)



# apply genotyping error and subsample to half the individuals
TW <- tweak_slim_simmed_data(
  X = simmed,
  S = 0.5, # use a small fraction for current tests
  EV = 0.01,
  ED = 0.004,
  MV = 0.25,
  MD = 0.25
)

# identify and extend ancestral segments
AS <- ancestral_n_segments(
  D = TW$spp_diag_snps,
  chrom_levels = unique(TW$spp_diag_snps$chrom)
)


## this block is removed and we use the MUCH faster bedtools version, now,
## which also keeps it all single threaded, which will be way better for
## efficient parallelization of the simulations
# future::plan(future::multicore, workers = 8)  # must set according to SLURM cores
#system.time(
#  E <- extend_ancestral_segments_3(
#    D = AS,
#    diag_spp_levels = c(tweaked$spp_names, "YCT")
#  )
#)


# get the ancestral segments using bedtools
E <- extend_ancestral_segments_using_bedtools(
  AS,
  species_levels = c(TW$spp_names, "YCT")
)


# Critical: We have to change the names of the species to A, B, C
# in the diagnostic species data frame
new_diag <- TW$spp_diag_snps %>%
  mutate(
    tmpf = as.integer(factor(diag_spp, levels = c(TW$spp_names, "YCT"))),
    diag_spp = c("A", "B", "C")[tmpf]
  ) %>%
  select(diag_spp, chrom, pos, indiv, n)

# Prepare the data structures
PRP <- prepare_for_gpp(
  V = TW$variable_snps,
  M = new_diag, #TW$spp_diag_snps,
  E = E,
  diag_err_rate = 0.005
)


# look at the ped_ids in different components of the simmed indivs
simmed$trueQ_and_pedigree %>% group_by(ind_pop, ind_time) %>% summarise(min = min(ped_id), max = max(ped_id))

# from that, we can get kids in different parts of the sim
inds_pop1_time0 <- {x <- unique(TW$variable_snps$indiv); x[as.integer(x) >= 20000 & as.integer(x) <= 20999]}
inds_pop1_time1 <- {x <- unique(TW$variable_snps$indiv); x[as.integer(x) >= 18000 & as.integer(x) <= 18999]}
inds_pop2_time0 <- {x <- unique(TW$variable_snps$indiv); x[as.integer(x) >= 21000 & as.integer(x) <= 21999]}
inds_pop2_time1 <- {x <- unique(TW$variable_snps$indiv); x[as.integer(x) >= 19000 & as.integer(x) <= 19999]}

#  compare pop1 kids from time 0 to parents at time 1
boink <- parallel::mclapply(c(inds_pop1_time0, inds_pop2_time0), function(kid) {
  print(kid)
  pairwise_genotype_probs(
    L = PRP,
    kid = kid,
    par = c(inds_pop1_time1, inds_pop2_time1)
  )
}, mc.cores = 8) %>%
  bind_rows()

got_em  <- boink %>%
  mutate(
    logl_ratio = probKidParental - probKidUnrel,
    kid_id_int = as.integer(kid_id)
  ) %>%
  arrange(desc(logl_ratio)) %>%
  left_join(simmed$trueQ_and_pedigree %>% filter(anc_pop == 1), by = join_by(kid_id_int == ped_id)) %>%
  mutate(is_a_parent = as.integer(par_id) == ped_p1 | as.integer(par_id) == ped_p2)


ggplot(got_em, aes(x = logl_ratio, fill = is_a_parent)) +
  geom_histogram(binwidth = 1) +
  xlim(-10, 125)



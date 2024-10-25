# deal with the logging
if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log$log, open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

library(tidyverse)
library(MixedUpParents)


# get the parameter values, and list some default ones for testing
# interactively (not within snakemake)
if(exists("snakemake")) {
  inrds  <- snakemake@input$inrds
  outrds <- snakemake@output$outrds
  threads <- snakemake@threads[[1]]
} else {
  inrds  <- "results/scenario-nonWF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/tweaked2mup.rds"
  outfile <- "test-mup-logls.rds"
  threads <- 8
}


TW <- read_rds(inrds)


# identify and extend ancestral segments
AS <- ancestral_n_segments(
  D = TW$spp_diag_snps,
  chrom_levels = unique(TW$spp_diag_snps$chrom)
)


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


# get the ids of the final generation samples and compare those
# to everyone else (including the non-self, final-generation samples)

# this preps a data frame
qped <- TW$sampled_pedQ %>%
  filter(anc_pop == 1)  # get just one row per indiv

# we are going to compare individuals to candidate parents in both populations,
#pops_and_years <- qped %>%
#  left_join(qped %>% select(ped_id, ind_pop) %>% rename(ped_p1_pop = ind_pop), by = join_by(ped_p1 == ped_id)) %>%
#  left_join(qped %>% select(ped_id, ind_pop) %>% rename(ped_p2_pop = ind_pop), by = join_by(ped_p2 == ped_id))

by_pop_all <- qped %>%
  group_by(ind_pop) %>%
  summarise(all_ids = list(unique(ped_id)))

last_generation_by_pop <- qped %>%
  filter(ind_time == 0) %>%
  group_by(ind_pop) %>%
  summarise(last_gen_ids = list(unique(ped_id)))

all_them <- last_generation_by_pop %>%
  left_join(by_pop_all, by = join_by(ind_pop))


# here is a function to run MUP when you pass it a vector of offspring and a vector
# of parental candidates
run_mup <- function(kids, pars) {
  parallel::mclapply(kids, function(kid) {
    candi <- setdiff(pars, kid)
    pairwise_genotype_probs(
      L = PRP,
      kid = kid,
      par = candi
    )
  }, mc.cores = threads) %>%
    bind_rows()
}



logls <- all_them %>%
  mutate(mup_output = map2(.x = last_gen_ids, .y = all_ids, .f = run_mup))

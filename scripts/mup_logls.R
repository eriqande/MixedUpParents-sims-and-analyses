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
  outrds <- "test-mup-logls.rds"
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

kids <- qped %>%
  filter(ind_time == 0) %>%
  pull(ped_id) %>%
  unique()

all_candi <- qped %>%
  pull(ped_id) %>%
  unique()


# now, we also want to keep a record of whether the parents were in the sample
# or not.
pars_in_samples <- qped %>%
  filter(ped_id %in% kids) %>%
  select(ped_id, ped_p1, ped_p2) %>%
  mutate(
    p1_in_sample = ped_p1 %in% all_candi,
    p2_in_sample = ped_p2 %in% all_candi,
    num_parents_in_sample = p1_in_sample + p2_in_sample
  ) %>%
  mutate(
    ped_id = as.character(ped_id),
    ped_p1 = as.character(ped_p1),
    ped_p2 = as.character(ped_p2)
  )

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



logls <- run_mup(kids, all_candi) %>%
  mutate(logl_ratio = probKidParental - probKidUnrel) %>% arrange(desc(logl_ratio))


# now, attach the age and pop information to each of these pair members
meta <- qped %>%
  select(ped_id, ind_pop, ind_time, admix_fract) %>%
  rename(id = ped_id, pop = ind_pop, time = ind_time, q1 = admix_fract) %>%
  mutate(id = as.character(id))

kid_meta <- meta
names(kid_meta) <- str_c("kid", names(meta), sep = "_")

par_meta <- meta
names(par_meta) <- str_c("par", names(meta), sep = "_")

# join that stuff on there, and also join on the relationship information.
# note that takes some ordering of the id's to make sure we have it right.
all_pairs <- logls %>%
  left_join(kid_meta, by = join_by(kid_id)) %>%
  left_join(par_meta, by = join_by(par_id)) %>%
  mutate(
    min_id = as.character(pmin(as.integer(kid_id), as.integer(par_id))),
    max_id = as.character(pmax(as.integer(kid_id), as.integer(par_id)))
  ) %>%
  left_join(
    TW$pairwise_relats %>%
      select(id_1, id_2, dom_relat, max_hit),
    by = join_by(min_id == id_1, max_id == id_2)
  ) %>%
  select(-min_id, -max_id) %>%
  mutate(
    dom_relat = replace_na(dom_relat, "U"),
    max_hit = replace_na(max_hit, 0L)
  )


# that is all the pairs, but, I don't want to write them all out.
# It will be more than sufficient to write out the top 10 parents for
# _from_each_cohort_ for each kid, and we will add in the information
# about whether the parents were sampled
trimmed_pairs <- all_pairs %>%
  arrange(kid_id, par_time, desc(logl_ratio))
  group_by(kid_id, par_time) %>%
  slice(1:10) %>%
  ungroup() %>%
  arrange(kid_id, logl_ratio) %>%
  left_join(
    pars_in_samples,
    by = join_by(kid_id == ped_id),
    relationship = "many-to-one"
  )



# that is just what we need to make the ROC curves.  So, write it out.
write_rds(trimmed_pairs, file = outrds, compress = "xz")

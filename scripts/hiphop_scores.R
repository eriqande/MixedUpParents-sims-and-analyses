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
  marker_set <- snakemake@params$marker_set  # could be var or both
} else {
  inrds  <- "test-tweaked2mup.rds"
  outrds <- "test-hot-scores.rds"
  threads <- 1
  marker_set <- "both"
}



IN <- read_rds(inrds)



# get the marker set to work with.  The user/workflow chooses whether to use
# just the variable SNPs or to use all the diagnostic SNPs in addition to the
# variable SNPs.
if(marker_set == "both") {
  MS <- bind_rows(
    IN$variable_snps,
    IN$spp_diag_snps %>% select(chrom, pos, indiv, n)
  ) %>%
    arrange(indiv, chrom, pos)
} else if(marker_set == "var") {
  MS <- IN$variable_snps
} else {
  stop("Unknown marker_set parameter: ", marker_set)
}


# now, we turn that dude into a matrix with positions along the rows and
# samples along the columns
genos_wide <- MS %>%
  pivot_wider(names_from = indiv, values_from = n)



# get the 0-based index of all the markers in a tibble
marker_indexes <- genos_wide %>%
  select(chrom, pos) %>%
  mutate(mIdx = 0:(n() - 1), .before = chrom)


# convert genos_wide into an integer matrix with just the markers
# and convert the NAs to -1
geno_mat <- genos_wide %>%
  select(-chrom, -pos) %>%
  as.matrix()
geno_mat[ is.na(geno_mat) ] <- -1L


# then get the 0-based indexes of the samples into a tibble
sample_indexes <- tibble(
  ped_id = colnames(geno_mat)
) %>%
  mutate(
    sIdx = 0:(n() - 1),
    .before = ped_id
  )



# this preps a data frame with the pedigree and admixture fractions of the individuals
qped <- IN$sampled_pedQ %>%
  filter(anc_pop == 1)  # get just one row per indiv


# we are going to compare individuals to candidate parents in both populations
# in this case we are going to want to get the 0-based indexes of the kids and then
# of all the candidates, but we also want vectors of their ped_ids
kids_tib <- qped %>%
  filter(ind_time == 0) %>%
  mutate(ped_id = as.character(ped_id)) %>%
  left_join(sample_indexes, by = join_by(ped_id))

kid_idxs <- kids_tib %>%
  pull(sIdx) %>%
  unique()

# this vector stores the ped_ids as characters
kids <- kids_tib %>%
  pull(ped_id) %>%
  unique()


all_candi_tib <- qped %>%
  mutate(ped_id = as.character(ped_id)) %>%
  left_join(sample_indexes, by = join_by(ped_id))

all_candi_idxs <- all_candi_tib %>%
  pull(sIdx) %>%
  unique()

# and here is the vector of ped_ids
all_candi <- all_candi_tib %>%
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



# here is a function to get the HOT statistics when you pass it a vector of offspring indexes, kids, and a vector
# of parental candidate indexes, pars
run_HOT <- function(kids, pars) {
  parallel::mclapply(kids, function(kid) {
    candi <- setdiff(pars, kid)
    calc_HOT(
      G = geno_mat,
      kid = kid,
      par = candi
    ) %>%
      as_tibble()
  }, mc.cores = threads) %>%
    bind_rows()
}



# now we should be able to run that
HOTs <- run_HOT(kid_idxs, all_candi_idxs) %>%
  mutate(hot_fract = num_hot / num_non_missing) %>%
  arrange(kIdx, hot_fract)



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
all_pairs <- HOTs %>%
  left_join(sample_indexes %>% rename(kid_id = ped_id), by = join_by(kIdx == sIdx)) %>%
  left_join(sample_indexes %>% rename(par_id = ped_id), by = join_by(pIdx == sIdx)) %>%
  left_join(kid_meta, by = join_by(kid_id)) %>%
  left_join(par_meta, by = join_by(par_id)) %>%
  mutate(
    min_id = as.character(pmin(as.integer(kid_id), as.integer(par_id))),
    max_id = as.character(pmax(as.integer(kid_id), as.integer(par_id)))
  ) %>%
  left_join(
    IN$pairwise_relats %>%
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
  arrange(kid_id, par_time, hot_fract) %>%
  group_by(kid_id, par_time) %>%
  slice(1:10) %>%
  ungroup() %>%
  arrange(kid_id, hot_fract) %>%
  left_join(
    pars_in_samples,
    by = join_by(kid_id == ped_id),
    relationship = "many-to-one"
  )

write_rds(trimmed_pairs, file = outrds, compress = "xz")



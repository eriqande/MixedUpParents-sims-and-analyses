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
  inrds  <- "results/scenario-nonWF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/tweaked2mup.rds"
  outrds <- "test-mup-logls.rds"
  threads <- 8
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
    )
  }, mc.cores = threads) %>%
    bind_rows()
}



# deal with the logging
if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log$log, open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

library(tidyverse)
library(MixedUpSlimSims)


# get the parameter values, and list some default ones for testing
# interactively (not within snakemake)
if(exists("snakemake")) {
  inrds  <- snakemake@input$inrds
  ppn_sampled <- as.numeric(snakemake@params$ppn_sampled)
  var_err <- as.numeric(snakemake@params$var_err)
  diag_err <- as.numeric(snakemake@params$diag_err)
  var_miss <- as.numeric(snakemake@params$var_miss)
  diag_miss <- as.numeric(snakemake@params$diag_miss)
  outfile <- snakemake@output$outrds
} else {
  inrds  <- "test-slim-sim.rds"
  ppn_sampled <- 0.5
  var_err <- 0.01
  diag_err <- 0.004
  var_miss <- 0.25
  diag_miss <- 0.25
  outfile <- "test-tweak2mup.rds"
}


simmed <- read_rds(inrds)

# apply genotyping error and subsampling and missing data
TW <- tweak_slim_simmed_data(
  X = simmed,
  S = ppn_sampled, # use a small fraction for current tests
  EV = var_err,
  ED = diag_err,
  MV = var_miss,
  MD = diag_miss
)


# now, we are going to add the pedigree and admixture fractions for the
# ped_ids that were actually sampled (note that some of the ped_p1 and ped_p2
# indivs in this tibble will not have been sampled...)
sample_ids <- unique(TW$variable_snps$indiv)
samp_ped <- simmed$trueQ_and_pedigree %>%
  filter(ped_id %in% sample_ids)

TW$sampled_pedQ <- samp_ped


write_rds(TW, file = outfile, compress = "xz")

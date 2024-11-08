# deal with the logging
if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log$log, open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

library(tidyverse)


# get the parameter values, and list some default ones for testing
# interactively (not within snakemake)
if(exists("snakemake")) {
  inList  <- snakemake@input$inList
  outrds <- snakemake@output$outrds
} else {
  inList <- list(
    "results/scenario-WF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/mup_rocs.rds",
    "results/scenario-WF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-1/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/mup_rocs.rds",
    "results/scenario-WF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-2/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/mup_rocs.rds"
  )
  outrds <- "results/scenario-nonWF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/mup_rocs.rds"
}


# pretty straightforward.  Just read the files in, put the path on there,
# bind_rows() them all, then parse the paths.
big_tib <- lapply(inList, function(x) read_rds(x) %>% mutate(path = x) %>% select(path, everything())) %>%
  bind_rows()

big_tib %>%
  extract(
    path,
    into = c(
      "scenario",
      "popsize1",
      "popsize2",
      "mig_rate1",  # fraction of pop2 that ends up migrating into pop1
      "mig_rate2",
      "rep",
      "samp_ppn",
      "var_err_rate",
      "diag_err_rate",
      "var_miss_rate",
      "diag_miss_rate",
      "remainder"
    ),
    regex = "results/scenario-(.+)/ps1-(.+)-ps2-(.+)-mr1-(.+)-mr2-(.+)/rep-(.+)/ppn-(.+)-verr-(.+)-derr-(.+)-vmiss-(.+)-dmiss-(.+)/(.+)$",
    remove = FALSE,
    convert = TRUE
  )  %>%
  write_rds(file = outrds, compress = "xz")

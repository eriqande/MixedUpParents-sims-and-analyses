




# this is in development

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
  inrds  <- snakemake@input$inrds
  outrds <- snakemake@output$outrds
} else {
  inrds  <- "results/scenario-nonWF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/mup_all_pairs.rds"
  outfile <- "results/scenario-nonWF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/mup_rocs.rds"
}


trimmed_pairs <- read_rds(inrds)

# after filtering out different candidate parent sets, this will compute the ROC
make_roc_mup <- function(X) {
  top_2s <- X %>%
    arrange(desc(logl_ratio)) %>%
    group_by(kid_id) %>%
    slice(1:2) %>%
    ungroup() %>%
    arrange(desc(logl_ratio))

  # now, we count the max number of correct parentage assignments
  tot_trues <- top_2s %>%
    distinct(kid_id, num_parents_in_sample) %>%
    pull(num_parents_in_sample) %>%
    sum()

  # and calculate the fraction of true-positives achieved, as well as the
  # fraction of false positives
  roc_values <- top_2s %>%
    mutate(
      tpr = cumsum(dom_relat == "PO") / tot_trues,
      num_false = cumsum(dom_relat != "PO"),
      fpr = num_false / max(num_false)
    )

  roc_values
}


# now we can do the two groups we want to
ret <- list(
  exclude_same_cohort = trimmed_pairs %>%
    filter(par_time > 0) %>%
    make_roc_mup(),
  include_same_cohort = trimmed_pairs %>%
    make_roc_mup()
) %>%
  bind_rows(.id = "cohort_inclusion")



write_rds(ret, file = outrds)

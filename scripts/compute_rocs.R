# Compute the ROC curves from the MUP and HipHop outputs

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
  #inrds  <- "results/scenario-nonWF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/mup_all_pairs.rds"
  inrds <- "results/scenario-nonWF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/hot_only_var_all_pairs.rds"
  outrds <- "results/scenario-nonWF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/mup_rocs.rds"
}


trimmed_pairs <- read_rds(inrds)

# after filtering out different candidate parent sets, this will compute the ROC
# this is set up so that it uses a metric that increases for pairs that are
# more likely (or have fewer mendelian incompatibilities).  We also
# set up what kind of inference procedure we used.
make_roc_mup <- function(X) {
  if("hot_fract" %in% names(X)) {
    X2 <- X %>%
      mutate(
        method = "hot",
        metric = -hot_fract,
        .before = kIdx
      )

  } else if("logl_ratio" %in% names(X)) {
    X2 <- X %>%
      mutate(
        method = "mup",
        metric = logl_ratio,
        .before = kIdx
      )

  } else {
    stop("Didn't find hot_fract or logl_ratio amongst the names")
  }

  top_2s <- X2 %>%
    arrange(desc(metric)) %>%
    group_by(kid_id, metric) %>%  # The lines up to the unnest are some fancy stuff to deal
                                  # with metrics that are the same.  If the metric is the same
                                  # we just randomly sort the individuals.
    nest() %>%
    mutate(
      data2 = map(.x = data, .f = function(x) x[sample(1:nrow(x)), ]),
      data = data2
    ) %>%
    select(-data2) %>%
    unnest(cols = data) %>%
    group_by(kid_id) %>%
    slice(1:2) %>%
    ungroup() %>%
    arrange(desc(metric))

  # now, we count the max number of correct parentage assignments
  tot_trues <- top_2s %>%
    distinct(kid_id, num_parents_in_sample) %>%
    pull(num_parents_in_sample) %>%
    sum()

  # get the total number of candidate offspring, too
  tot_candi_kid <- n_distinct(top_2s$kid_id)

  # and calculate the fraction of true-positives achieved, as well as the
  # fraction of false positives
  roc_values <- top_2s %>%
    mutate(
      tpr = cumsum(dom_relat == "PO") / tot_trues,
      num_false = cumsum(dom_relat != "PO"),
      fpr = num_false / (tot_candi_kid * 2.0)  # each kid has two parents to assign
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

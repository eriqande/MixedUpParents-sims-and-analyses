# deal with the logging
if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log$log, open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
}


library(tidyverse)
library(sequoia)



# get the parameter values, and list some default ones for testing
# interactively (not within snakemake)
if(exists("snakemake")) {
  inrds <- snakemake@input$inrds
  seq_sim_results <- snakemake@output$seq_sim_results
  seq_roc_results <- snakemake@output$seq_roc_results
  which_markers <- snakemake@params$which_markers
  cohort_situation <- snakemake@params$cohort_situation
} else {
  inrds <- "development/small-langford-tweaked2mup.rds"
  seq_sim_results <- "sequoia_simres.rds"
  seq_roc_results <- "sequoia_rocres.rds"
  which_markers <- "only_var"  # can be either "only_var" or "both_diag_and_var"
  cohort_situation <- "exclude_same_cohort" # can be either exclude_same_cohort or include_same_cohort
}

##################################
##### Read in the input files #####
##################################
ind_info <- read_rds(inrds)

########################################
##### Put data into Sequoia format #####
########################################

if(which_markers == "only_var") {
  MarkerTib <- ind_info$variable_snps
} else if(which_markers == "both_diag_and_var") {
  MarkerTib <- ind_info$spp_diag_snps %>%
    select(-diag_spp) %>%
    bind_rows(ind_info$variable_snps)
} else {
  stop("Unrecognized which_markers: ", which_markers, ".  Must be either only_var or both_diag_and_var.")
}

diag_and_var_snps.tib <- MarkerTib %>%
  mutate(across(everything(), ~replace(., is.na(.), "-9"))) %>%
  pivot_wider(names_from = c(chrom, pos), values_from = n) %>%
  mutate_all(as.numeric)

# A function to turn a tibble into a matrix
make_matrix <- function(tib, rownames = NULL){
  my_matrix <-  as.matrix(tib)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}

seq.mat <- make_matrix(select(diag_and_var_snps.tib,-indiv),
                       pull(diag_and_var_snps.tib, indiv))




##### Make the life history data frame #####

seq_metadata.df <- ind_info$sampled_pedQ %>% select(ped_id, ind_time, sex) %>%
  distinct() %>% #Get rid of one duplicate row per individual
  mutate(sex = replace(sex, sex == 0, 2)) %>% #Change sex code from 0 to 2
  mutate(ind_time = case_when(ind_time == 0 ~ 11, #Change birth years
                              ind_time == 1 ~ 10,
                              ind_time == 2 ~ 9)) %>%
  rename(ID = ped_id, BirthYear = ind_time, Sex = sex) %>% #Rename the columns
  as.data.frame()


# some notes here:
# the ped_p1 parent is given as Sex = 2, which means male in Sequoia
# the ped_p2 parent is given as Sex = 1, which means female in Sequoia



####################################
##### Run the Sequoia analysis #####
####################################

# Currently it is always exclude_same_cohort.  We need to add a tweak to
# do include_same_cohort when that is what we are doing...

sequoia_par_run <- sequoia(GenoM = seq.mat, LifeHistData = seq_metadata.df,
                           MaxSibshipSize = 500,
                           Module = "par", Plot = F, Err = 1e-02,
                           Tassign = 0,
                           Tfilter = -2e100,
                           quiet = "verbose")


write_rds(sequoia_par_run, file = seq_sim_results, compress = "xz")

############################################
##### Compare sequoia results to truth
#############################################

# Eric threw this in here to be closer to how it was done for MUP and hiphop,
# and to just go ahead and write out the ROC curve for it.

# First thing we do is get the sequoia results listed on a per-parent basis.
# If a pair of parents was assigned, we given both parents the same LLRPair
# log-likelihood ratio.
seq_ped_results2 <- sequoia_par_run$PedigreePar %>%
  as_tibble() %>%
  rename(
    kid_id = id,
    LLR_dam = LLRdam,
    LLR_sire = LLRsire
  ) %>%
  mutate(
    LLR_dam = ifelse(!is.na(LLRpair), LLRpair, LLR_dam),
    LLR_sire = ifelse(!is.na(LLRpair), LLRpair, LLR_sire)
  ) %>%
  select(-OHdam, -OHsire, -MEpair)

seq_ped_results3 <- seq_ped_results2 %>%
  pivot_longer(c(dam, sire), names_to = "par_type", values_to = "par_id") %>%
  mutate(logl_ratio = ifelse(par_type == "dam", LLR_dam, LLR_sire)) %>%
  select(kid_id, par_id, par_type, logl_ratio)


# now, we just need to join stuff on there to get the kid and parent
# pops and times and q1's.  We get a tibble of that extra information
basic_meta <- ind_info$sampled_pedQ %>%
  filter(anc_pop == 1) %>%
  select(ped_id, ind_pop, ind_time, admix_fract) %>%
  rename(ind_q1 = admix_fract) %>%
  mutate(ped_id = as.character(ped_id))

# and we will also get the true pedigree in terms of dams and sires
dam_sire_ped <- ind_info$sampled_pedQ %>%
  select(starts_with("ped")) %>%
  distinct() %>%
  rename(ped_sire = ped_p1, ped_dam = ped_p2) %>%
  pivot_longer(c(ped_sire, ped_dam)) %>%
  mutate(name = str_replace(name, "ped_", "")) %>%
  rename(par_type = name, true_parent = value) %>%
  mutate(ped_id = as.character(ped_id))

# and here is a vector of all the sampled individuals
sampled_indivs <- ind_info$sampled_pedQ %>%
  pull(ped_id) %>%
  unique() %>%
  as.character()

ready_for_roc <- seq_ped_results3 %>%
  left_join(basic_meta, by = join_by(kid_id == ped_id)) %>%
  rename(kid_pop = ind_pop, kid_time = ind_time, kid_q1 = ind_q1)  %>%
  left_join(basic_meta, by = join_by(par_id == ped_id)) %>%
  rename(par_pop = ind_pop, par_time = ind_time, par_q1 = ind_q1) %>%
  left_join(dam_sire_ped, by = join_by(kid_id ==ped_id, par_type)) %>%
  mutate(
    parent_correct = par_id == true_parent,
    parent_in_sample = true_parent %in% sampled_indivs
  ) %>%
  arrange(desc(logl_ratio)) %>%
  mutate(cohort_inclusion = cohort_situation, .before = kid_id)


# and now we compute the ROC parts.  Note that we are giving the
# false positive rate as a fraction of the max number of false
# positives.  This might break if there are none because sequoia
# is too conservative about assigning anyone.

# Get the total number of True Positives possible
tot_poss_tp <- sum(ready_for_roc$parent_in_sample)

# And also the total number of parents that could be assigned
# to the candidate kids
tot_candi_kid_par <- nrow(ready_for_roc)


# We are going to let the FPR be the fraction of all candidate kids made up to a
# certain point (likelihood) that are incorrect.
roc <- ready_for_roc %>%
  mutate(
    tpr = cumsum(parent_correct) / tot_poss_tp,
    num_false = cumsum(parent_correct == FALSE),
    fpr = num_false / tot_candi_kid_par
  )


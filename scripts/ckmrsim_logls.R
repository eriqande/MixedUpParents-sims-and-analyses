# this script just treats the admixed samples as all coming from a single
# population and calculates pairwise logls for parentage and full sibs
# from that, and then chooses the best. This is probably a better comparison
# than Sequoia, which is incredibly slow.


# deal with the logging
if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log$log, open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

library(tidyverse)
library(CKMRsim)

# get the parameter values, and list some default ones for testing
# interactively (not within snakemake)
if(exists("snakemake")) {
  inrds  <- snakemake@input$inrds
  outrds <- snakemake@output$outrds
  outnumLoc <- snakemake@output$outnumLoc
  marker_set <- snakemake@params$marker_set
  threads <- snakemake@threads[[1]]
} else {
  inrds  <- "results/scenario-nonWF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/tweaked2mup.rds"
  outrds <- "test-ckmr-logls.rds"
  marker_set <- "only_var"
  threads <- 8
}


TW <- read_rds(inrds)


# Make a long-genos data frame suitable for CKMRsim from the
# variable and also the diagnostic SNPs. Let's do that with a function.
make_long_genos <- function(X) {
  # get rid of the diag_spp column if it is there
  X2 <- X[, names(X) != "diag_spp"]

  X2 %>%
    rename(
      Indiv = indiv,
      Chrom = chrom,
      Pos = pos
    ) %>%
    mutate(
      Locus = str_c(Chrom, Pos, sep = ":"),
      alle_1 = case_match(n, 0 ~ "A", 1 ~ "A", 2 ~ "B"),
      alle_2 = case_match(n, 0 ~ "A", 1 ~ "B", 2 ~ "B")
    ) %>%
    select(-n) %>%
    pivot_longer(cols = c(alle_1, alle_2), values_to = "Allele", names_to = c("out", "gene_copy"), names_sep = "_") %>%
    select(-out, -Chrom, -Pos)
}


var_long_genos <- make_long_genos(TW$variable_snps)
diag_long_genos <- make_long_genos(TW$spp_diag_snps)

if(marker_set == "both_diag_and_var") {
  long_genos <- bind_rows(var_long_genos, diag_long_genos)

} else if(marker_set == "only_var") {
  long_genos <- var_long_genos
} else {
  stop("Unkown value of marker_set: ", marker_set)
}



# Now, get the allele freqs from the long_genos. Some chicanery here because reindex markers doesn't work with multiple real chromosomes
locus_names <- unique(long_genos$Locus)
afreqs_ready <- long_genos %>%
  count(Locus, Allele) %>%
  group_by(Locus) %>%
  mutate(
    Freq = n / sum(n),
    Chrom = "Unk",
    Pos = as.integer(factor(Locus, levels = locus_names))
  ) %>%
  ungroup() %>%
  select(Chrom, Pos, Locus, Allele, Freq) %>%
  arrange(Chrom, Pos, desc(Freq)) %>%
  mutate(AlleIdx = NA,
         LocIdx = NA) %>%
  filter(!is.na(Allele)) %>%
  reindex_markers()



# make a CKMR object
ckmr <- create_ckmr(
  D = afreqs_ready,
  kappa_matrix = kappas[c("PO", "FS", "U"), ],
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = 0.005),
  ge_mod_true_pars_list = list(epsilon = 0.005)
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


# Calculate all the logl ratios we need
pw_2_LRTs <- lapply(
  X = list(
    POU = c("PO", "U"),
    POFS = c("PO", "FS")
  ),
  FUN = function(x) {
    pairwise_kin_logl_ratios(
      D1 = long_genos %>% filter(Indiv %in% all_candi),
      D2 = long_genos %>% filter(Indiv %in% kids),
      CK = ckmr,
      numer = x[1],
      denom = x[2],
      num_cores = threads,
    )
  }
) %>%
  bind_rows(
    .id = "lr_type"
  ) %>%
  pivot_wider(names_from = lr_type, values_from = logl_ratio) %>%
  arrange(D2_indiv, desc(POU)) %>%
  filter(D2_indiv != D1_indiv)  # we have to remove comparisons to themselves


# now, we prepare this to look more like something that might have come out of MUP
logls <- pw_2_LRTs %>%
  mutate(
    kid_id = D2_indiv,
    par_id = D1_indiv,
    logl_ratio = POU
  )



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
  arrange(kid_id, par_time, desc(logl_ratio)) %>%
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

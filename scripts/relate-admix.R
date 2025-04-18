# this script gets the simulated genotypes, then it
# puts them in plink binary format and runs ADMIXTURE on them
# and then it runs them through RelateAdmix.


# deal with the logging
if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log$log, open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

library(tidyverse)
library(vroom)

# get the parameter values, and list some default ones for testing
# interactively (not within snakemake)
if(exists("snakemake")) {
  inrds  <- snakemake@input$inrds
  outrds <- snakemake@output$outrds
  plinkped <- snakemake@output$plinkped
  plinkbed <- snakemake@output$plinkbed
  indsfile <- snakemake@output$indsfile
  output_k <- snakemake@output$output_k  # this is what relateAdmix spits out
  threads <- snakemake@threads[[1]]
} else {
  inrds  <- "results/scenario-nonWF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/tweaked2mup.rds"
  outrds <- "test-ckmr-logls.rds"
  threads <- 8
  plinkped <- "plink.ped"
  plinkbed <- "binary.bed"
  indsfile <- "indsfile.tsv"
  output_k <- "output.k"
}


TW <- read_rds(inrds)


# put all the genotypes (diag and var) together, and make the genotypes
# two columns, like "A\tG" for a het (1), "A\tA" for 0, and "G\tG" for 2.
genos <- bind_rows(
  TW$variable_snps,
  TW$spp_diag_snps %>% select(-diag_spp)
) %>%
  arrange(indiv, chrom, pos) %>%
  mutate(
    geno2 = case_match(
      n,
      0L ~ "A\tA",
      1L ~ "A\tG",
      2L ~ "G\tG"
    )
  ) %>%
  mutate(
    geno2 = replace_na(geno2, "0\t0")
  )

# get the chroms and pos
map <- genos %>%
  distinct(chrom, pos) %>%
  mutate(
    id = sprintf("Locus%04d", 1:n()),
    junk = 0L,
    .after = chrom
  ) %>%
  mutate(
    chrom = as.integer(str_replace_all(chrom, "[^0-9]", ""))
  )

wide <- pivot_wider(
  genos %>% select(-n),
  names_from = c(chrom, pos),
  values_from = geno2,
  names_sep = ":"
) %>%
  mutate(
    indiv = str_c("x", indiv),
  ) %>%
  mutate(
    d1 = 0,
    d2 = 0,
    d3 = 0,
    d4 = -9,
    .after = indiv
  ) %>%
  mutate(
    pop = "pop1",
    .before = indiv
  )



# write that out to the plink file
write.table(
  wide,
  file = plinkped,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  map,
  file = plinkmap,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)


# also, write out the inds file:
wide %>%
  select(indiv) %>%
  mutate(
    indiv = str_replace(indiv, "x", ""),
    idx = (1:n()) - 1L
  ) %>%
  write_tsv(., file = indfile)



# now, we want to make a binary fileset out of that
CALL <- paste(
  "plink --file",
  str_replace(plinkped, "\\.ped$", ""),
  "--aec --make-bed --biallelic-only --chr-set 29 --out",
  str_replace(plinkbed, "\\.bed$", "")
)

system(CALL)

# The following two steps will be done in a temp directory
# because both ADMIXTURE and relateAdmix seem to just write
# their output files out in the current directory.  We will then
# grab the results from those and copy them over or read them
# in and process them as necessary.


# now, we need to run admixture on this
CALL2 <- paste(
  "admixture",
  plinkbed,
  "2"
)
system(CALL2)


# relateAdmix -plink binary.bed -f binary.2.P -q binary.2.Q -P 8

#### We also want to get some backstory on all the samples  ####
# this stuff is analogous/identical to what gets done in mup_logls.R
# this preps a data frame
qped <- TW$sampled_pedQ %>%
  filter(anc_pop == 1)  # get just one row per indiv

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




#### Read in the results ####

# We want to read in the results and have them be very similar in structure
# to the mup_logls results so that we can analyse them in the same way.

# relateAdmix compares pairs only once with low index in column one and
# higher index in column 2.  mup_logls, on the other hand, compares each
# kid in the first (kid_id) column to all others.  So, we are going to have to
# modify the output of relateAdmix so as to:
#  1. add rows with the pairs listed in reverse order
#  2. retain only those rows that have a kid listed in the first column
#
# throughout that we also want to keep a min_id and a max_id column to be sure
# we are joining correctly onto the pairwise relationships tibble.

# then read in the results and join the ids back on there, and sort
# the id columns so the lower value is always first. (it already is, but
# just in case...) and add the true relationshps on there
inds_tib <- read_tsv(indsfile)
ra_output <- vroom(output_k) %>%
  left_join(inds_tib %>% rename(indiv1 = indiv), by = join_by(ind1 == idx) ) %>%
  left_join(inds_tib %>% rename(indiv2 = indiv), by = join_by(ind2 == idx) ) %>%
  mutate(
    min_id = as.character(pmin(as.integer(indiv1), as.integer(indiv2))),
    max_id = as.character(pmax(as.integer(indiv1), as.integer(indiv2)))
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


# just gonna look at these values for a moment

# stick the pairwise relationships on there:


g <- ggplot(
  ra_output %>% mutate(relat = str_c(dom_relat, max_hit, sep = "-")),
  aes(x = k1, y = k2, colour = relat)) +
  geom_point() +
  facet_wrap(~relat)


ggsave(g, filename = "scatter.jpg", width = 14, height = 14)





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

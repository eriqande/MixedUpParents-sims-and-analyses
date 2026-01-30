# this script gets the simulated genotypes, then it
# puts them in plink binary format and runs ADMIXTURE on them
# and then it runs them through RelateAdmix.

#saveRDS(snakemake, file = "/tmp/snakemake.rds")
#stop()


# deal with the logging
if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log$Rlog, open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

library(tidyverse)
library(vroom)

# get the parameter values, and list some default ones for testing
# interactively (not within snakemake)
if(exists("snakemake")) {
  inrds  <- snakemake@input$inrds
  marker_set <- snakemake@params$marker_set
  outrds <- snakemake@output$outrds
  plinkped <- snakemake@output$plinkped
  plinkbed <- snakemake@output$plinkbed
  plinkmap <- snakemake@output$plinkmap
  indsfile <- snakemake@output$indsfile
  admixture_plinkped <- snakemake@output$admixture_plinkped
  admixture_plinkbed <- snakemake@output$admixture_plinkbed
  admixture_plinkmap <- snakemake@output$admixture_plinkmap
  relate_admix_locus_mask <- snakemake@output$relate_admix_locus_mask
  indsfile <- snakemake@output$indsfile
  output_k <- snakemake@output$output_k  # this is what relateAdmix spits out
  binaryP <- snakemake@output$binaryP  # this is what relateAdmix uses from ADMIXTURE
  binaryQ <- snakemake@output$binaryQ  # this is what relateAdmix uses from ADMIXTURE
  Rlog <- snakemake@log$Rlog
  admixture_log <- snakemake@log$admixture_log
  relateAdmix_log <- snakemake@log$relateAdmix_log
  plink_log <- snakemake@log$plink_log
  admixture_plink_log <- snakemake@log$admixture_plink_log
  munge_call_log <- snakemake@log$munge_call_log
  threads <- snakemake@threads[[1]]
  marker_set <- snakemake@params$marker_set
} else {
  inrds  <- "results/scenario-nonWF_simple/ps1-1200-ps2-1200-mr1-0.06-mr2-0.02/rep-0/ppn-0.5-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/tweaked2mup.rds"
  marker_set <- "only_var"
  outrds <- "test-ckmr-logls.rds"
  threads <- 8
  plinkped <- "plink.ped"
  plinkbed <- "binary.bed"
  plinkmap <- "plink.map"
  admixture_plinkped <- "admixture_plink.ped"
  admixture_plinkbed <- "admixture_binary.bed"
  admixture_plinkmap <- "admixture_plink.map"
  relate_admix_locus_mask <- "relate_admix_locus_mask.txt"
  indsfile <- "indsfile.tsv"
  output_k <- "output.k"
  Rlog <- "Rlog"
  plink_log <- "plink_log"
  admixture_plink_log <- "admixture_plink_log"
  admixture_log <- "admixture_log"
  relateAdmix_log <- "relateAdmix_log"
  munge_call_log <- "munge_call_log"
}


TW <- read_rds(inrds)


############ here, we have to choose which marker_set we are using (only_var vs both_diag_and_var) ###########
#### Actually, it is a little more complicated.  We should use the diagnostic markers (or maybe all the markers) in order
#### to run admixture on them, then we can explore whether we want the only_vars or both the variable
#### and diagnostic ones to be used for the relateAdmix.  Doing so is going to required parsing the ADMIXTURE .P results
#### and keeping only the variable markers, when relateAdmix is not using the diagnostic markers
#### for itself (i.e. when doing only_var).

# put all the genotypes (diag and var) together, and make the genotypes
# two columns, like "A\tG" for a het (1), "A\tA" for 0, and "G\tG" for 2.
if(marker_set == "only_var") {
  genos_input <- TW$variable_snps
} else if(marker_set == "both_diag_and_var") {
  genos_input <- bind_rows(
    TW$variable_snps,
    TW$spp_diag_snps %>% select(-diag_spp)
  )
} else {
  stop("Unrecognized marker set: ", marker_set)
}

########
# here we are getting all the markers in the admixture data set
admixture_genos_input <- bind_rows(
  TW$variable_snps,
  TW$spp_diag_snps %>% select(-diag_spp)
)

admixture_genos <- admixture_genos_input %>%
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
admixture_map <- admixture_genos %>%
  distinct(chrom, pos) %>%
  mutate(
    id = sprintf("Locus%04d", 1:n()),
    junk = 0L,
    .after = chrom
  ) %>%
  mutate(
    chrom = as.integer(str_replace_all(chrom, "[^0-9]", ""))
  )

admixture_wide <- pivot_wider(
  admixture_genos %>% select(-n),
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

# Working here to output an ADMIXTURE-specific file with both marker sets
# write that out to the plink file
write.table(
  admixture_wide,
  file = admixture_plinkped,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  admixture_map,
  file = admixture_plinkmap,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

#######

# Here we are getting the genotypes that relateAdmix will use
genos <- genos_input %>%
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



# write that out to the plink file for relateAdmix use
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
  write_tsv(., file = indsfile)



# now, we want to make a binary fileset out of both the Admixture and the RelateAdmix
# versions
CALL <- paste(
  "eval ' plink --file",
  str_replace(plinkped, "\\.ped$", ""),
  "--aec --make-bed --biallelic-only --chr-set 29 --out",
  str_replace(plinkbed, "\\.bed$", ""),
  " > ", plink_log, " 2>&1 ' "
)

ADMIXTURE_CALL <- paste(
  "eval ' plink --file",
  str_replace(admixture_plinkped, "\\.ped$", ""),
  "--aec --make-bed --biallelic-only --chr-set 29 --out",
  str_replace(admixture_plinkbed, "\\.bed$", ""),
  " > ", admixture_plink_log, " 2>&1 ' "
)

system(CALL)
system(ADMIXTURE_CALL)


# Now, we make a 0/1 mask for which lines in the allele frequency file that ADMIXTURE
# outputs should be retained for use in relateAdmix.  If relateAdmix is using both markers
# sets (var and diagnostic) then all the markers should be used.  If relateAdmix is
# using just the variable markers, then only the those are used.  The way we make this
# 0/1 mask is just by recording which loci from admixture_wide are found in the relateAdmix wide genos.
admixture_locus_names <- names(admixture_wide)[-(1:6)]
relateadmix_locus_names <- names(wide)[-(1:6)]

ra_locus_01_mask <- admixture_locus_names %in% relateadmix_locus_names %>%
  as.integer()

cat(ra_locus_01_mask, sep = "\n", file = relate_admix_locus_mask)

# The following two steps will be done in a temp directory
# because both ADMIXTURE and relateAdmix seem to just write
# their output files out in the current directory.  We will then
# grab the results from those and copy them over or read them
# in and process them as necessary.


tmpdir <- tempfile()
dir.create(tmpdir, recursive = TRUE)
CurrentWD <- getwd()  # use for making absolute paths


message("Running admixture and relateAdmix in ", tmpdir)

# note that all paths will be relative to the top-level project directory,
# so we can  recreate that as needed
admixture_plinkbed_absolute <- file.path(CurrentWD, admixture_plinkbed)
plinkbed_absolute <- file.path(CurrentWD, plinkbed)
admixture_log_absolute <- file.path(CurrentWD, admixture_log)
relateAdmix_log_absolute <- file.path(CurrentWD, relateAdmix_log)
relate_admix_locus_mask_absolute <- file.path(CurrentWD, relate_admix_locus_mask)
munge_call_log_absolute <- file.path(CurrentWD, munge_call_log)
# now, we need to run admixture on this
ADMIXTURE_CALL2 <- paste(
  " cd ", tmpdir,
  "; eval ' admixture ",
  admixture_plinkbed_absolute,
  " 2 -j", threads, " > ", admixture_log_absolute, " 2>&1 ' ",
  sep = ""
)

# this runs admixture in the tmpdir and the results "binary.2.P" and "binary.2.Q"
# get written in that tmpdir.
system(ADMIXTURE_CALL2)


# now, we create binary.2.P and binary.2.Q as needed for relateAdmix to work
MUNGE_CALL <- paste(
  " cd ", tmpdir, "; cp admixture_binary.2.Q binary.2.Q",
  "; eval ' paste ", relate_admix_locus_mask_absolute,
  "  admixture_binary.2.P | grep ^1 | cut -f 2,3 > binary.2.P 2> ", munge_call_log_absolute, "  ' ",
  sep = ""
)


system(MUNGE_CALL)


# So, now, we run relateAdmix the same way.  This is the second place the multiple threads come in.
# It will get done fast
CALL3 <-  paste(
  " cd ", tmpdir,
  "; eval ' relateAdmix -plink",
  plinkbed_absolute,
  " -f binary.2.P -q binary.2.Q -P ",
  threads,
  " > ", relateAdmix_log_absolute, " 2>&1 ' "
)

system(CALL3)


# Now, we need to copy output.k to the correct location
stopifnot(file.copy(from = file.path(tmpdir, basename(output_k)), to = output_k, overwrite = TRUE) == TRUE)

# and we also need to copy the binaryP and binaryQ files that were used as input to
# RelateAdmix (they were output from ADMIXTURE, but also the .P file had rows
# removed if it is only_var).
stopifnot(file.copy(from = file.path(tmpdir, basename(binaryP)), to = binaryP, overwrite = TRUE) == TRUE)
stopifnot(file.copy(from = file.path(tmpdir, basename(binaryQ)), to = binaryQ, overwrite = TRUE) == TRUE)


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



#### Read in the results ####

# We want to read in the results and have them be very similar in structure
# to the mup_logls results so that we can analyse them in the same way.

# relateAdmix compares pairs only once with low index in column one and
# higher index in column 2.  mup_logls, on the other hand, compares each
# kid in the first (kid_id) column to all others.  So, we are going to have to
# modify the output of relateAdmix so as to:
#  1. add rows with the pairs listed in reverse order
#  2. retain only those rows that have a kid listed in the first column


# then read in the results and join the ids back on there, and sort
# the id columns so the lower value is always first. (it already is, but
# just in case...) and add the true relationshps on there
inds_tib <- read_tsv(indsfile)
ra_output <- vroom(output_k) %>%
  left_join(inds_tib %>% rename(indiv1 = indiv), by = join_by(ind1 == idx) ) %>%
  left_join(inds_tib %>% rename(indiv2 = indiv), by = join_by(ind2 == idx) )



# now, we prepare this to look more like something that might have come out of MUP.
# Recall that we need to orient each pair so that the kid is first.  We do that by
# making an additional tibble that has all the IDs reversed, and then retaining only those
# that have kids in the first ID column.
ra_output_rev <- ra_output %>%
  mutate(
    tmp_ind = ind1,
    ind1 = ind2,
    ind2 = tmp_ind,
    tmp_indiv1 = indiv1,
    indiv1 = indiv2,
    indiv2 = tmp_indiv1
  )


ra_output2 <- bind_rows(ra_output, ra_output_rev) %>%
  filter(indiv1 %in% kids) %>%
  select(-starts_with("tmp"), -ind1, -ind2) %>%
  select(indiv1, indiv2, everything()) %>%
  mutate(
    indiv1 = as.character(indiv1),
    indiv2 = as.character(indiv2)
  )

# now, ra_output2 has all the stuff that we need. (I've checked the number of
# rows and it checks out:  equal to (num_kids * num_candi) - num_kids
# ra_output2 essentially serves the role of the tibble "logls" in the MUP calculations.
# Note that "indiv1" is the kid of each pair.

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
# note that takes some ordering of the id's to make sure we have it right:
# we must keep a min_id and a max_id column to be sure
# we are joining correctly onto the pairwise relationships tibble.
# And we also calculate a coefficient of relatedness, which is what we use
# as our metric for ordering pairs.  This is 0.5*k1 + k2.
all_pairs <- ra_output2 %>%
  left_join(kid_meta, by = join_by(indiv1 == kid_id)) %>%
  left_join(par_meta, by = join_by(indiv2 == par_id)) %>%
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
  ) %>%
  mutate(
    relate_admix_metric = k1,
    .after = indiv2
  )





# that is all the pairs, but, I don't want to write them all out.
# It will be more than sufficient to write out the top 10 parents for
# _from_each_cohort_ for each kid, and we will add in the information
# about whether the parents were sampled
trimmed_pairs <- all_pairs %>%
  arrange(indiv1, par_time, desc(relate_admix_metric)) %>%
  group_by(indiv1, par_time) %>%
  slice(1:10) %>%   # take the top 10 from each possible level of parent time
  ungroup() %>%
  left_join(
    pars_in_samples,
    by = join_by(indiv1 == ped_id),
    relationship = "many-to-one"
  ) %>%
  arrange(indiv1, desc(relate_admix_metric))




# that is just what we need to make the ROC curves.  So, write it out.
# Though right now I am writing out all pairs because I am trying to
# figure out what to use as a metric...
write_rds(trimmed_pairs, file = outrds, compress = "xz")

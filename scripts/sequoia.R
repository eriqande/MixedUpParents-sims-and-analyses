# deal with the logging
if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log$log, open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

library(sequoia)

# get the parameter values, and list some default ones for testing
# interactively (not within snakemake)
if(exists("snakemake")) {
  inrds <- snakemake@input$inrds
  slim_inrds <- snakemake@input$slim_inrds
  seq_results <- snakemake@output$seq_results
} else {
  inrds <- "~/Documents/UM/Research_Documents/Hybrid_Simulations/MUP-Sequoia/results/scenario-nonWF_simple/ps1-1000-ps2-1000-mr1-0.0-mr2-0.025/rep-0/ppn-0.25-verr-0.01-derr-0.004-vmiss-0.25-dmiss-0.25/tweaked2mup.rds"
  slim_inrds <- "~/Documents/UM/Research_Documents/Hybrid_Simulations/MUP-Sequoia/results/scenario-nonWF_simple/ps1-1000-ps2-1000-mr1-0.0-mr2-0.025/rep-0/slim-output.rds"
  seq_results <- "~/Documents/UM/Research_Documents/Hybrid_Simulations/MUP-Sequoia/results/seq-results.rds"
}

##################################
##### Read in the input files #####
##################################
ind_info <- read_rds(inrds)
slim_sim <- read_rds(slim_inrds)

########################################
##### Put data into Sequoia format #####
########################################

diag_and_var_snps.tib <- ind_info$spp_diag_snps %>% select(-diag_spp) %>%
  bind_rows(ind_info$variable_snps) %>%
  mutate(across(everything(), ~replace(., is.na(.), "-9"))) %>%
  pivot_wider(names_from = c(chrom, pos), values_from = n) %>%
  mutate_all(as.numeric)

#A function to turn a tibble into a matrix
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


####################################
##### Run the Sequoia analysis #####
####################################

sequoia_par_run <- sequoia(GenoM = seq.mat, LifeHistData = seq_metadata.df,
                           MaxSibshipSize = 500,
                           Module = "par", Plot = F,
                           quiet = "verbose")

#####################################################
##### Check the accuracy of the Sequoia results #####
#####################################################

#Get the true trio information
true_trios.tib <- ind_info$sampled_pedQ %>% filter(ind_time == 0, anc_pop == 1) %>% #Select the offspring generation, and only chromosome info. from pop 1 (RBT pop)
  select(ped_id, admix_fract, ped_p1, ped_p2, sex) %>%
  mutate(sex = replace(sex, sex == 0, 2)) #Change sex code from 0 to 2; 1 = F, 2 = M
  
### Assign the parents from SLiM to their respective sexes given in Sequoia ###
# Needed for parsing Sequoia results that are by dam and sire
# 1 = F, 2 = M
dam_list <- c()
sire_list <- c()

for (ind_row in 1:nrow(true_trios.tib)) {
  
  #Get IDs of the two parents
  p1_id <- true_trios.tib$ped_p1[ind_row]
  p2_id <- true_trios.tib$ped_p2[ind_row]
  
  #Get sex info. of parent 1 and add the individual either to dam or sire list
  p1_info <- slim_sim$trueQ_and_pedigree[which(slim_sim$trueQ_and_pedigree$ped_id == p1_id), ]
  
  if (p1_info$sex[1] == 1) {
    dam_list <- c(dam_list, p1_id)
  }
  
  if (p1_info$sex[1] == 0) {
    sire_list <- c(sire_list, p1_id)
  }
  
  #Get sex info. of parent 2 and add the individual either to dam or sire list
  p2_info <- slim_sim$trueQ_and_pedigree[which(slim_sim$trueQ_and_pedigree$ped_id == p2_id), ]
  
  if (p2_info$sex[1] == 1) {
    dam_list <- c(dam_list, p2_id)
  }
  
  if (p2_info$sex[1] == 0) {
    sire_list <- c(sire_list, p2_id)
  }
  
}

#Assign the dams and sires to their respective columns
#Remove the ped_p1 and ped_p2 columns and make a new data frame with dames and sires
true_trios_sexes.tib <- true_trios.tib %>%
  mutate(True.dam = dam_list, True.sire = sire_list) %>%
  select(-ped_p1, -ped_p2) %>%
  mutate(sex = case_when(sex == 1 ~ "F",
                         sex == 2 ~ "M"))


### Now check to see if Sequoia results were correct or not, LLR scores, and how many OH mismatches there were ###

#Add in columns that will be filled in the for-loop below
true_trios_sexes.tib$Inferred.dam <- NA
true_trios_sexes.tib$Inferred.sire <- NA
true_trios_sexes.tib$Dam.correct <- NA
true_trios_sexes.tib$Sire.correct <- NA
true_trios_sexes.tib$LLRdam <- NA
true_trios_sexes.tib$LLRsire <- NA
true_trios_sexes.tib$OHdam <- NA
true_trios_sexes.tib$OHsire <- NA
true_trios_sexes.tib$Dam.in.analysis <- NA #Record whether the dam was analyzed (sampled) or not
true_trios_sexes.tib$Sire.in.analysis <- NA #Record whether the sire was analyzed (sampled) or not

##### Get IDs of the individuals that were removed during subsampling #####
#All individuals in the analysis (after subsampling)
sample_ids <- unique(ind_info$variable_snps$indiv)

#All individuals in the pedigree (before subsampling) are in slim_sim$trueQ_and_pedigree
removed_inds.tib <- slim_sim$trueQ_and_pedigree %>%
  filter(ind_time <= 2, !(ped_id %in% sample_ids), anc_pop == 1)
  
removed_inds <- removed_inds.tib$ind_id


## Check the results to see how many offspring were correctly assigned
for (offsp_result in 1:nrow(true_trios_sexes.tib)) {
  
  #Get the true dam and sire for the offspring
  true_triad <- true_trios_sexes.tib[offsp_result, ]
  true_offsp_id <- true_triad$ped_id
  true_dam <- true_triad$True.dam
  true_sire <- true_triad$True.sire
  
  #Find whether the parents were analyzed/sampled or not
  ifelse(true_dam %in% removed_inds,
         true_trios_sexes.tib$Dam.in.analysis[offsp_result] <- 0,
         true_trios_sexes.tib$Dam.in.analysis[offsp_result] <- 1)
  
  ifelse(true_sire %in% removed_inds,
         true_trios_sexes.tib$Sire.in.analysis[offsp_result] <- 0,
         true_trios_sexes.tib$Sire.in.analysis[offsp_result] <- 1)
  
  #Get offspring, dam and sire ID from Sequoia results
  sim_result_triad <- sequoia_par_run$PedigreePar[which(sequoia_par_run$PedigreePar$id == true_offsp_id), ]
  
  sim_result_dam <- sim_result_triad$dam
  sim_result_sire <- sim_result_triad$sire
  
  #Write inferred dam and sire to the results data frame
  true_trios_sexes.tib$Inferred.dam[offsp_result] <- sim_result_dam
  true_trios_sexes.tib$Inferred.sire[offsp_result] <- sim_result_sire
  
  #Enter here if a dam has been assigned
  if (!is.na(sim_result_dam)) {
    
    if (true_dam == sim_result_dam) {
      
      true_trios_sexes.tib$Dam.correct[offsp_result] <- 1
      true_trios_sexes.tib$LLRdam[offsp_result] <- sim_result_triad$LLRdam
      true_trios_sexes.tib$OHdam[offsp_result] <- sim_result_triad$OHdam
      
    } else if (true_dam != sim_result_dam) { #Incorrect dam assignment
      
      true_trios_sexes.tib$Dam.correct[offsp_result] <- 0
      true_trios_sexes.tib$LLRdam[offsp_result] <- sim_result_triad$LLRdam
      true_trios_sexes.tib$OHdam[offsp_result] <- sim_result_triad$OHdam
      
    }
    
  }
  
  #Enter here if a sire has been assigned
  if (!is.na(sim_result_sire)) {
    
    if (true_sire == sim_result_sire) {
      
      true_trios_sexes.tib$Sire.correct[offsp_result] <- 1
      true_trios_sexes.tib$LLRsire[offsp_result] <- sim_result_triad$LLRsire
      true_trios_sexes.tib$OHsire[offsp_result] <- sim_result_triad$OHsire
      
    } else if (true_sire != sim_result_sire) { #Incorrect sire assignment
      
      true_trios_sexes.tib$Sire.correct[offsp_result] <- 0
      true_trios_sexes.tib$LLRsire[offsp_result] <- sim_result_triad$LLRsire
      true_trios_sexes.tib$OHsire[offsp_result] <- sim_result_triad$OHsire
      
    }
    
  }
  
  if (is.na(sim_result_dam)) { #Dames weren't assigned
    
    true_trios_sexes.tib$Dam.correct[offsp_result] <- 0
    true_trios_sexes.tib$LLRdam[offsp_result] <- NA
    true_trios_sexes.tib$OHdam[offsp_result] <- NA
    
  }
  
  if (is.na(sim_result_sire)) { #Sires weren't assigned
    
    true_trios_sexes.tib$Sire.correct[offsp_result] <- 0
    true_trios_sexes.tib$LLRsire[offsp_result] <- NA
    true_trios_sexes.tib$OHsire[offsp_result] <- NA
    
  }
  
}

write_rds(true_trios_sexes.tib, file = seq_results, compress = "xz")

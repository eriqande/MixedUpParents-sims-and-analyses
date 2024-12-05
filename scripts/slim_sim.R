
# deal with the logging
if(exists("snakemake")) {
  # redirect output and messages/errors to the log
  log <- file(snakemake@log$log, open="wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

library(tidyverse)
library(MixedUpSlimSims)
library(reticulate)


# get the parameter values, and list some default ones for testing
# interactively (not within snakemake)
if(exists("snakemake")) {
  slimfile  <- snakemake@input$slimfile
  varfreq_file  <- snakemake@input$varfreqs
  diagmarker_file  <- snakemake@input$diagmarkers
  pop_size1 <- as.integer(snakemake@params$pop_size1)
  pop_size2 <- as.integer(snakemake@params$pop_size1)
  mig_rate1 <- as.numeric(snakemake@params$mig_rate1)
  mig_rate2 <- as.numeric(snakemake@params$mig_rate2)
  outfile <- snakemake@output$outrds
} else {
  slimfile <- "config/slim_templates/WF_eric_test.slim"
  varfreq_file <- "config/var_freqs/cyclone-var-freqs.rds"
  diagmarker_file <- "config/diag_markers/wct-rbt-yct-spp-diag-markers.rds"
  pop_size1 <- 990L
  pop_size2 <- 980L
  mig_rate1 <- 0.055
  mig_rate2 <- 0.011
  outfile <- "test-slim-sim.rds"
}


# read the slim template in and replace the values as needed
slim_text <- readChar(slimfile, file.info(slimfile)$size) %>%
  str_glue(.open = "<<<<", .close = ">>>>")

# write that out to a temp file
slim_tmp <- tempfile()
cat(slim_text, file = slim_tmp)

# get the varfreqs and the diagmarkers
varfreqs <- read_rds(varfreq_file)
diagmarkers <- read_rds(diagmarker_file)

# get the seed for slim by hashing the output file name into an 8 digit number
seed <- rlang::hash(outfile) %>%
  str_replace_all("[^0-9]", "") %>%
  str_sub(end = 8) %>%
  as.integer()

# divide that seed by 2 to get the R-seed
Rseed <- floor(seed / 2)
set.seed(Rseed)

# do this explicitly so reticulate does not ask if you want to
# make an environment
use_condaenv(Sys.getenv("MUP_CONDA"))



# simulate data using MixedUpSlimSims
simmed <- slim_sim_a_dataset(
  AF = varfreqs,
  genome_info = mykiss_chroms,
  founder_pop_nums = c(pop_size1, pop_size2),
  slim_file = slim_tmp,
  slim_seed = seed,
  years_list = list(`1` = 0:10, `2` = 0:10),
  marker_pop_time_list = list(p1 = 9:11, p2 = 9:11),
  diagnostic_markers = diagmarkers
)

simmed$slim_seed <- seed
simmed$R_seed <- Rseed

write_rds(simmed, file = outfile, compress = "xz")


quit(save="no", status=0, runLast=FALSE)


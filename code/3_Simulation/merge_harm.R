rm(list = ls())

INTERACTIVE <- FALSE

library(optparse)
library(tidyverse)
library(devtools)

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="working directory"),
    make_option(c("-d", "--combat_dir"), action="store", default=NULL,
        type="character", help="output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (INTERACTIVE) {
    opt$combat_dir <- "simulation_revisions/other_models/sim_50_0.1_0.1_0.0"
}

setwd(opt$wdir)
load_all("submodules/byyfunctions")

harm_files <- list.files(path = opt$combat_dir, pattern = "sim_harm.*\\.csv", full.names = TRUE)
sim_harm <- harm_files %>% map(read_csv) %>% bind_rows()

write_csv(sim_harm, file.path(opt$combat_dir, "sim_harm.csv"))

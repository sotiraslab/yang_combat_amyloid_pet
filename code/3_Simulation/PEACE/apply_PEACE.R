# ==============================================================================

# Name: apply_PEACE.R
# Author: Braden Yang
# Created: 08/24/2024
# Description: Apply PEACE to harmonize simulated data

# ==============================================================================

rm(list = ls())

INTERACTIVE <- FALSE

library(optparse)
library(tidyverse)
library(rstan)
library(devtools)

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="working directory"),
    make_option(c("-o", "--odir"), action="store", default=NULL,
        type="character", help="output directory"),
    make_option(c("-p", "--peace_dir"), action="store", default="data/simulation/PEACE",
        type="character", help="directory containing trained PEACE models for simulation experiment"),
    make_option(c("-m", "--mode"), action="store", default=NULL,
        type="character", help="PEACE mode ('nocovar', 'age', 'age_sex_apoe')"),
    make_option(c("-f", "--filepath"), action="store", default=NULL,
        type="character", help="path to simulated data CSV")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (INTERACTIVE) {
    opt$mode <- "age_sex_apoe"
    opt$filepath <- "data/simulation/sim_10_0.1_0.1_0.0/sim.csv"
}

setwd(opt$wdir)
load_all("submodules/byyfunctions")
source("code/1_Harmonize/PEACE/peace_functions.R")

stan_dir <- "code/1_Harmonize/PEACE"
byyfunctions::make_dir(opt$odir)

# load simulated data
sim_df <- read_csv(opt$filepath)

# load PEACE model
if (opt$mode == "nocovar") {
    model_path <- file.path(opt$peace_dir, "model__nocovar.RDS")
} else if (opt$mode == "age") {
    model_path <- file.path(opt$peace_dir, "model__age.RDS")
} else {
    model_path <- file.path(opt$peace_dir, "model__age_sex_apoe.RDS")
}
model <- readRDS(model_path)
feature_names <- read_csv("data/csv/feature_names_trunc.csv", col_names = FALSE, show_col_types = FALSE) %>% pull(1)

# apply PEACE model
sim_data <- get_PEACE_data(sim_df, opt$mode)
sim_peace <- apply_PEACE_wrapper(model, sim_data, sim_df, str_glue("peace__{opt$mode}"), opt$mode == "nocovar")

# save results
write_csv(sim_peace, file.path(opt$odir, str_glue("sim_harm_PEACE__{opt$mode}.csv")))

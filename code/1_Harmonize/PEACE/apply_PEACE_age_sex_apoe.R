# ==============================================================================

# Name: apply_PEACE_age_sex_apoe.R
# Author: Braden Yang
# Created: 09/19/2024
# Description: apply trained PEACE model to head-to-head crossover data, age +
#   sex + APOE as covariates

# ==============================================================================

# clear environment
rm(list = ls())

library(tidyverse)
library(rstan)
library(optparse)

library(devtools)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="Path to project directory (if none, uses current working directory)"),
    make_option(c("-o", "--odir"), action="store", default="data/head2head/PEACE",
        type="character", help="Path to output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# set working directory
if (!is.null(opt$wdir)) {
    setwd(opt$wdir)
}

# load byyfunctions
devtools::load_all("submodules/byyfunctions")

# make odir
byyfunctions::make_dir(opt$odir)

# ============================
# ===== SOURCE FUNCTIONS =====
# ============================

source("code/1_Harmonize/PEACE/peace_functions.R")

# ==============================
# ===== LOAD AND PREP DATA =====
# ==============================

# load data
crossover_df <- read_csv("data/demographics/crossover.csv")
crossover_train_df <- read_csv("data/demographics/crossover_train.csv")
feature_names <- read_csv("data/csv/feature_names_trunc.csv", col_names = FALSE) %>% pull(1)

# extract train data
data_train__age_sex_apoe_outfile <- file.path(opt$odir, "data_train__age_sex_apoe.RDS")
if (file.exists(data_train__age_sex_apoe_outfile)) {
    data_train__age_sex_apoe <- read_rds(data_train__age_sex_apoe_outfile)
} else {
    data_train__age_sex_apoe <- get_PEACE_data(crossover_train_df, "age_sex_apoe")
    write_rds(data_train__age_sex_apoe, data_train__age_sex_apoe_outfile)
}

# extract test (crossover) data
data_crossover__age_sex_apoe_outfile <- file.path(opt$odir, "data_crossover__age_sex_apoe.RDS")
if (file.exists(data_crossover__age_sex_apoe_outfile)) {
    data_crossover__age_sex_apoe <- read_rds(data_crossover__age_sex_apoe_outfile)
} else {
    data_crossover__age_sex_apoe <- get_PEACE_data(crossover_df, "age_sex_apoe")
    write_rds(data_crossover__age_sex_apoe, data_crossover__age_sex_apoe_outfile)
}

# =====================
# ===== FIT PEACE =====
# =====================

model__age_sex_apoe_outfile <- file.path(opt$odir, "model__age_sex_apoe.RDS")
if (file.exists(model__age_sex_apoe_outfile)) {
    model__age_sex_apoe <- read_rds(model__age_sex_apoe_outfile)
} else {
    model__age_sex_apoe <- stan(
        file = "code/1_Harmonize/PEACE/PEACE_mult_covar.stan",
        data = data_train__age_sex_apoe,
        seed = 42,
        iter = 100
    )
    model__age_sex_apoe@stanmodel@dso <- new("cxxdso")
    write_rds(model__age_sex_apoe, model__age_sex_apoe_outfile)
}

# =======================
# ===== APPLY PEACE =====
# =======================

crossover_peace__age_sex_apoe <- apply_PEACE_wrapper(model__age_sex_apoe, data_crossover__age_sex_apoe, crossover_df, "peace__age_sex_apoe")
write_rds(crossover_peace__age_sex_apoe, file.path(opt$odir, "crossover_PEACE__age_sex_apoe.RDS"))

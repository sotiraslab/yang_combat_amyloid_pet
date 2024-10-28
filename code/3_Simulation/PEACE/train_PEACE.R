# ==============================================================================

# Name: train_PEACE.R
# Author: Braden Yang
# Created: 08/24/2024
# Description: Train PEACE model to harmonize simulated data

# ==============================================================================

rm(list = ls())

library(optparse)
library(tidyverse)
library(rstan)
library(devtools)

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="working directory"),
    make_option(c("-o", "--odir"), action="store", default=NULL,
        type="character", help="output directory"),
    make_option(c("-m", "--mode"), action="store", default=NULL,
        type="character", help="PEACE mode ('nocovar', 'age', 'age_sex_apoe')")
)
opt <- parse_args(OptionParser(option_list = option_list))

setwd(opt$wdir)
load_all("submodules/byyfunctions")
source("code/1_Harmonize/PEACE/peace_functions.R")

stan_dir <- "code/1_Harmonize/PEACE"
byyfunctions::make_dir(opt$odir)

# load training data
train_df <- read_csv("data/simulation/amypos.csv")
feature_names <- read_csv("data/csv/feature_names_trunc.csv", col_name = FALSE) %>% pull(1)

if (opt$mode == "nocovar") {
    # ===== PEACE - NO COVARIATES =====
    data_outfile <- file.path(opt$odir, "data__nocovar.RDS")
    data_train <- get_PEACE_data(train_df, "nocovar")
    write_rds(data_train, data_outfile)
    
    model_outfile <- file.path(opt$odir, "model__nocovar.RDS")
    model <- stan(
        file = file.path(stan_dir, "PEACE_no_covar.stan"),
        data = data_train,
        seed = 42,
        iter = 100
    )
    model@stanmodel@dso <- new("cxxdso")
    write_rds(model, model_outfile)
} else if (opt$mode == "age") {
    # ===== PEACE - AGE =====
    data_outfile <- file.path(opt$odir, "data__age.RDS")
    data_train <- get_PEACE_data(train_df, "age")
    write_rds(data_train, data_outfile)

    model_outfile <- file.path(opt$odir, "model__age.RDS")
    model <- stan(
        file = file.path(stan_dir, "PEACE_mult_covar.stan"),
        data = data_train,
        seed = 42,
        iter = 100
    )
    model@stanmodel@dso <- new("cxxdso")
    write_rds(model, model_outfile)
} else {
    # ===== PEACE - AGE, SEX, APOE =====
    data_outfile <- file.path(opt$odir, "data__age_sex_apoe.RDS")
    data_train <- get_PEACE_data(train_df, "age_sex_apoe")
    write_rds(data_train, data_outfile)

    model_outfile <- file.path(opt$odir, "model__age_sex_apoe.RDS")
    model <- stan(
        file = file.path(stan_dir, "PEACE_mult_covar.stan"),
        data = data_train,
        seed = 42,
        iter = 100
    )
    model@stanmodel@dso <- new("cxxdso")
    write_rds(model, model_outfile)
}

# ==============================================================================

# Name: get_excel.R
# Author: Braden Yang
# Created: 02/14/2024
# Description: create excel spreadsheets with data

# ==============================================================================

# clear environment
rm(list = ls())

# ===========================
# ===== IMPORT PACKAGES =====
# ===========================

library(tidyverse)
library(gt)
library(optparse)

library(devtools)
library(openxlsx)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="Path to project directory (if none, uses current working directory)"),
    make_option(c("-o", "--odir"), action="store", default="tables",
        type="character", help="Path to output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# set working directory
if (!is.null(opt$wdir)) {
    setwd(opt$wdir)
}

# load byyfunctions
load_all("submodules/byyfunctions")

# create odir
byyfunctions::make_dir(opt$odir)

# =====================
# ===== LOAD DATA =====
# =====================

icc_df <- read_rds("data/head2head/icc.RDS")
mae_df <- read_rds("data/head2head/mae.RDS")
icc_ttest <- read_rds("data/head2head/icc_ttest.RDS")
mae_ttest <- read_rds("data/head2head/mae_ttest.RDS")
sim_df <- read_rds("data/simulation/stats.RDS")

# =========================
# ===== SAVE AS EXCEL =====
# =========================

l <- list(
    "ICC" = icc_df,
    "MAE" = mae_df,
    "ICC_ttest" = icc_ttest,
    "MAE_ttest" = mae_ttest,
    "simulation" = sim_df
)
openxlsx::write.xlsx(l, file = file.path(opt$odir, "stats.xlsx"))

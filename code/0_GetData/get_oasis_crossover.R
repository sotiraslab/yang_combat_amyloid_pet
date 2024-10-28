# ==============================================================================

# Name: get_crossover.R
# Author: Braden Yang
# Created: 03/13/2023
# Description: get crossover subjects from OASIS, who have both AV45 and PIB
# amyloid PET within a short period of time

# ==============================================================================

# clear environment
rm(list = ls())

# *** toggle variables ***
INTERACTIVE <- FALSE
# *** toggle variables ***

# ===========================
# ===== IMPORT PACKAGES =====
# ===========================

library(tidyverse)
library(optparse)
library(devtools)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="Path to project directory (if none, uses current working directory)"),
    make_option(c("-o", "--odir"), action="store", default="data",
        type="character", help="Path to output directory"),
    make_option(c("-p", "--pvc"), action="store_true", default=FALSE,
        type="logical", help="Extract partial volume corrected SUVRs from OASIS"),
    make_option(c("-m", "--white_matter"), action="store_true", default=FALSE,
        help="extract WM regions for ComBat harmonization")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (INTERACTIVE) {
    opt$odir <- "data"
    opt$white_matter <- TRUE
}

# set working directory
if (!is.null(opt$wdir)) {
    setwd(opt$wdir)
}

# load byyfunctions
load_all("submodules/byyfunctions")

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

remove_double_match <- function(.data, time_col1, time_col2, group_col) {

    .data <- .data %>%
        group_by({{group_col}}) %>%
        slice_min(abs({{time_col1}} - {{time_col2}})) %>%
        ungroup()

    return(.data)

}

select_feature <- function(av45_col, pib_col, tracer_col) {

    v <- case_when(
        {{tracer_col}} == "AV45" ~ {{av45_col}},
        {{tracer_col}} == "PIB"  ~ {{pib_col}},
        .default = NA
    )

    return(v)

}

# =====================
# ===== LOAD DATA =====
# =====================

if (opt$white_matter) {
    pet_df <- read_rds("data/pet_wm.RDS")
} else {
    pet_df <- read_rds("data/pet.RDS")
}
feature_names <- read_csv("data/csv/feature_names.csv", col_names = FALSE, show_col_types = FALSE) %>% pull(1)
if (opt$white_matter) {
    wm_feature_names <- read_csv("data/csv/white_matter_regions.csv", col_names = FALSE, show_col_types = FALSE) %>% pull(1)
    feature_names <- c(feature_names, wm_feature_names)
}
select_cols <- c(feature_names, "cdr", "mmse", "amyloid_positive", "clinical_group", "clinical_group_extended", "idx")

oasis_df <- pet_df %>%
    filter(
        study == "OASIS",
        pvc == opt$pvc
    )

# extract OASIS AV45 and PIB dataframes
oasis_av45_df <- oasis_df %>%
    filter(tracer == "AV45") %>%
    rename_with(~paste0(.x, ".AV45"), all_of(select_cols))
oasis_pib_df <- oasis_df %>%
    filter(tracer == "PIB") %>%
    rename_with(~paste0(.x, ".PIB"), all_of(select_cols))

# ===============================
# ===== GET CROSSOVER SCANS =====
# ===============================

# match every AV45 with PIB
crossover_df <- 
    byyfunctions::match_data(
        oasis_av45_df,
        match_df = oasis_pib_df,
        ref_col = "age",
        match_col = "age",
        select_col = paste0(select_cols, ".PIB"),
        group_col = "subj",
        max_diff = 90 / 365.25
    ) %>%
    filter(!is.na(age_match_age)) %>%  # select only successful matches
    rename(age.AV45 = age, age.PIB = age_match_age) %>%  # rename columns
    rename_with(~str_remove(.x, "age_match_"), starts_with("age_match_")) %>%
    select(-tracer) %>%  # drop previous tracer column
    mutate(pair_id = row_number()) %>%  # col to keep track of AV45-PIB pairing
    remove_double_match(age.AV45, age.PIB, idx.PIB)  # remove instances of AV45-to-PIB double-matching

# tidy data
crossover_tidy <- crossover_df %>%
    pivot_longer(
        cols = starts_with(feature_names),
        names_pattern = "(.*)\\.(AV45|PIB)",
        names_to = c("feature_roi", "tracer"),
        values_to = "feature_value"
    ) %>%
    pivot_wider(
        names_from = "feature_roi",
        values_from = "feature_value"
    ) %>%
    mutate(
        age = select_feature(age.AV45, age.PIB, tracer),
        cdr = select_feature(cdr.AV45, cdr.PIB, tracer),
        mmse = select_feature(mmse.AV45, mmse.PIB, tracer),
        amyloid_positive = select_feature(amyloid_positive.AV45, amyloid_positive.PIB, tracer),
        clinical_group = select_feature(clinical_group.AV45, clinical_group.PIB, tracer),
        clinical_group_extended = select_feature(clinical_group_extended.AV45, clinical_group_extended.PIB, tracer),
        idx = select_feature(idx.AV45, idx.PIB, tracer)
    )

# =====================
# ===== SAVE DATA =====
# =====================

# create output directory
byyfunctions::make_dir(opt$odir, verbose = FALSE)

if (opt$white_matter) {
    out_suffix <- "_wm"
} else {
    out_suffix <- ""
}

print("++ Saving RDS file ++")
filename <- if (opt$pvc) {str_glue("crossover_pvc{out_suffix}.RDS")} else {str_glue("crossover{out_suffix}.RDS")}
write_rds(crossover_tidy, file.path(opt$odir, filename))

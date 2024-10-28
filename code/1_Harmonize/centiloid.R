# ==============================================================================

# Name: centiloid.R
# Author: Braden Yang
# Created: 06/23/2023
# Description: compute centiloid (global and regional) using equations provided
#   by ADNI and OASIS datasets

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
        type="character", help="Path to output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# set working directory
if (!is.null(opt$wdir)) {
    setwd(opt$wdir)
}

# load byyfunctions
load_all("submodules/byyfunctions")

# make output directory
byyfunctions::make_dir(opt$odir)

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

adni_av45_to_centiloid <- function(av45_suvr) {

    # this equation was obtained from the Centiloid Level-2
    # analysis of the ADNI PET pipeline, written by the Pitt
    # group (although they used the processing of the UC
    # Berkeley group)

    cl <- 196.9 * av45_suvr - 196.03

    return(cl)

}

adni_fbb_to_centiloid <- function(fbb_suvr) {

    # this equation was obtained from the Centiloid Level-2
    # analysis of the ADNI PET pipeline, written by the Pitt
    # group (although they used the processing of the UC
    # Berkeley group)

    cl <- 159.08 * fbb_suvr - 151.65

    return(cl)

}

oasis_av45_to_centiloid <- function(av45_suvr, pvc = FALSE) {

    # these equations were obtained from the KARI methods and
    # definitions document for data freeze 15 (20190423)

    if (pvc) {
        cl <- 53.6 * av45_suvr - 43.2
    } else {
        cl <- 163.6 * av45_suvr - 181.0
    }
    return(cl)

}

oasis_pib_to_centiloid <- function(pib_suvr, pvc = FALSE) {

    # these equations were obtained from the KARI methods and
    # definitions document for data freeze 15 (20190423)

    if (pvc) {
        cl <- 45.0 * pib_suvr - 47.5
    } else {
        cl <- 111.8 * pib_suvr - 119.3
    }
    return(cl)

}

compute_centiloid <- function(suvr, tracer, study, pvc) {

    centiloid <- case_when(
        {{tracer}} == "AV45" & {{study}} == "ADNI" ~ adni_av45_to_centiloid({{suvr}}),
        {{tracer}} == "FBB" & {{study}} == "ADNI" ~ adni_fbb_to_centiloid({{suvr}}),
        {{tracer}} == "AV45" & {{study}} == "OASIS" & !pvc ~ oasis_av45_to_centiloid({{suvr}}),
        {{tracer}} == "PIB" & {{study}} == "OASIS" & !pvc ~ oasis_pib_to_centiloid({{suvr}}),
        {{tracer}} == "AV45" & {{study}} == "OASIS" & pvc ~ oasis_av45_to_centiloid({{suvr}}, pvc = TRUE),
        {{tracer}} == "PIB" & {{study}} == "OASIS" & pvc ~ oasis_pib_to_centiloid({{suvr}}, pvc = TRUE),
    )

    return(centiloid)

}

centiloid_across <- function(.data) {

    .data <- .data %>% 
        mutate(across(
            .cols = all_of(feature_names),
            .fns = ~ compute_centiloid(.x, tracer, study, pvc)
        ))

}

# =====================
# ===== LOAD DATA =====
# =====================

pet_df <- read_rds("data/pet.RDS")
baseline_df <- read_rds("data/baseline.RDS")
baseline_pvc_df <- read_rds("data/baseline_pvc.RDS")
crossover_df <- read_rds("data/crossover.RDS")
crossover_pvc_df <- read_rds("data/crossover_pvc.RDS")
feature_names <- read_csv("data/csv/feature_names.csv", col_names = FALSE) %>% pull(1)

# =============================
# ===== COMPUTE CENTILOID =====
# =============================

centiloid_list <- list()
centiloid_list[["centiloid"]] <- pet_df %>% centiloid_across()
centiloid_list[["centiloid_baseline"]] <- baseline_df %>% centiloid_across()
centiloid_list[["centiloid_baseline_pvc"]] <- baseline_pvc_df %>% centiloid_across()
centiloid_list[["centiloid_crossover"]] <- crossover_df %>% centiloid_across()
centiloid_list[["centiloid_crossover_pvc"]] <- crossover_pvc_df %>% centiloid_across()

# ===========================
# ===== SAVE CENTILOIDS =====
# ===========================

for (key in names(centiloid_list)) {
    write_rds(
        centiloid_list[[key]],
        file.path(opt$odir, paste0(key, ".RDS"))
    )
}

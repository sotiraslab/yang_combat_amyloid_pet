# ==============================================================================

# Name: apply_long_combat.R
# Author: Braden Yang
# Created: 08/24/2024
# Description: Train and apply longCombat (https://github.com/jcbeer/longCombat)
#   to simulated data

# ==============================================================================

rm(list = ls())

INTERACTIVE <- FALSE

library(optparse)
library(dplyr)
library(forcats)
library(readr)
library(purrr)
library(stringr)
library(longCombat)
library(devtools)

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="working directory"),
    make_option(c("-o", "--odir"), action="store", default=NULL,
        type="character", help="output directory"),
    make_option(c("-f", "--filepath"), action="store", default=NULL,
        type="character", help="path to simulated data CSV")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (INTERACTIVE) {
    opt$odir <- "data/revisions/simulation/longCombat"
    opt$mode <- "nocovar"
    opt$filepath <- "data/simulation/sim_10_0.1_0.1_0.0/sim.csv"
}

setwd(opt$wdir)
load_all("submodules/byyfunctions")

byyfunctions::make_dir(opt$odir)

harmonize_long_combat <- function(.data, group_id, formula_str) {
    # apply longCombat
    harm <- longCombat(
        data = .data,
        idvar = "subj",
        timevar = "age",
        batchvar = "tracer.fct",
        features = feature_names,
        formula = formula_str,
        ranef = "(1|subj)",
        verbose = FALSE
    )

    # return harmonized data
    harm_tibble <- harm$data_combat %>%
        as_tibble() %>%
        mutate(sim_idx = group_id[[1]])
    return(harm_tibble)
}

# load simulated data
sim_df <- read_csv(opt$filepath)
sim_train <- sim_df %>%
    mutate(
        tracer.fct = as_factor(tracer),
        sex_M = case_match(sex, "M"~1, "F"~0),
    )
feature_names <- read_csv("data/csv/feature_names_trunc.csv", col_names = FALSE, show_col_types = FALSE) %>% pull(1)

# apply longCombat to each simulation individually (age, sex, APOE)
sim_longCombat_temp <- sim_train %>%
    group_by(sim_idx) %>%
    group_map(~ harmonize_long_combat(.x, .y, "")) %>%
    reduce(bind_rows) %>%
    rename_with(~str_replace(.x, ".combat", ""), starts_with(feature_names))
sim_longCombat__nocovar <- dplyr::full_join(
    sim_df %>% select(-all_of(feature_names)),
    sim_longCombat_temp %>% select(sim_idx, subj, age, all_of(feature_names)),
    by = c("sim_idx", "subj", "age")
) %>%
    mutate(harmonization_method = "longCombat__nocovar")

# apply longCombat to each simulation individually (age, sex, APOE)
sim_longCombat_temp <- sim_train %>%
    group_by(sim_idx) %>%
    group_map(~ harmonize_long_combat(.x, .y, "age + sex_M + apoe")) %>%
    reduce(bind_rows) %>%
    rename_with(~str_replace(.x, ".combat", ""), starts_with(feature_names))
sim_longCombat__age_sex_apoe <- dplyr::full_join(
    sim_df %>% select(-all_of(feature_names)),
    sim_longCombat_temp %>% select(sim_idx, subj, age, all_of(feature_names)),
    by = c("sim_idx", "subj", "age")
) %>%
    mutate(harmonization_method = "longCombat__age_sex_apoe")

# save
write_csv(
    sim_longCombat__nocovar,
    file.path(opt$odir, "sim_harm_longCombat__nocovar.csv")
)
write_csv(
    sim_longCombat__age_sex_apoe,
    file.path(opt$odir, "sim_harm_longCombat__age_sex_apoe.csv")
)

# ==============================================================================

# Name: 4_statistical_test.R
# Author: Braden Yang
# Created: 12/15/2023
# Description: Perform statistical significance testing on simulated data

# ==============================================================================

rm(list = ls())

# *** toggle variables ***
INTERACTIVE <- FALSE
# *** toggle variables ***

# ===========================
# ===== IMPORT PACKAGES =====
# ===========================

library(optparse)
library(tidyverse)

library(lme4)
library(broom)
library(broom.mixed)

library(devtools)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="Path to project directory (if none, uses current working directory)"),
    make_option(c("-i", "--input"), action="store", default=NULL,
        type="character", help="path to sim_harm.csv"),
    make_option(c("-o", "--odir"), action="store", default="data/simulation",
        type="character", help="path to output directory"),
    make_option(c("-f", "--feature_roi"), action="store", type="character", default="summary",
        help="name of FS ROI to perform simulations on"),
    make_option(c("-t", "--tracer_covariate"), action="store_true", default=FALSE,
        help="include tracer as a covariate in LME model to test for significant differences"),
    make_option(c("-u", "--unharmonized"), action="store_true", default=FALSE,
        help="indicates that input is the unharmonized data ('sim.csv')")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (INTERACTIVE) {
    opt$feature_roi <- "summary"
    opt$tracer_covariate <- TRUE
}

# set working directory
if (!is.null(opt$wdir)) {
    setwd(opt$wdir)
}

# load byyfunctions
load_all("submodules/byyfunctions")

# make output directory
byyfunctions::make_dir(opt$odir)

# ============================
# ===== DEFINE VARIABLES =====
# ============================

# linear mixed effects formulas
resp_var <- "feature_value"
fixed_eff <- c("1", "age.baseline", "time.from.baseline", "sex", "apoe")
if (opt$tracer_covariate) {
    fixed_eff <- c(fixed_eff, "tracer")
}
treatment_eff <- c("trial.group", "time.from.baseline_trial.group")
random_eff <- "(1 | subj)"
lhs <- paste(c(fixed_eff, treatment_eff, random_eff), collapse = " + ")
lmer_formula <- paste(c(resp_var, "~", lhs), collapse = " ")
term_to_test <- "time.from.baseline_trial.group"

# get list of regions
roi_list <- read_csv("data/csv/feature_names_trunc.csv", col_names = FALSE) %>% pull(1)
roi_remove <- roi_list[!(roi_list %in% opt$feature_roi)]

# =====================
# ===== LOAD DATA =====
# =====================

sim_df <- read_csv(opt$input)
if (opt$unharmonized) {
    sim_df <- sim_df %>% mutate(harmonization_method = "suvr")
}
sim_df <- sim_df %>%
    select(-all_of(roi_remove)) %>%
    pivot_longer(  # pivot ROIs
        cols = all_of(opt$feature_roi),
        names_to = "feature_roi",
        values_to = "feature_value"
    )

# ====================
# ===== FIT LMER =====
# ====================

group_col <- c("sim_idx", "harmonization_method", "feature_roi")

lmer_df <- sim_df %>%
    nest(.by = all_of(group_col)) %>%
    mutate(
        lmer_model = map(
            data,
            ~ lmer(
                formula = formula(lmer_formula),
                data = .x
            )
        ),
        lmer_tidy = map(lmer_model, tidy),
        lrt = map(
            lmer_model,
            ~ drop1(.x, test = "Chisq") %>% tidy()
        )
    )

# extract model params and LRT results
param_df <- lmer_df %>%
    select(all_of(group_col), lmer_tidy) %>%
    unnest(lmer_tidy) %>%
    filter(term == term_to_test)
lrt_df <- lmer_df %>%
    select(all_of(group_col), lrt) %>%
    unnest(lrt) %>%
    filter(term == term_to_test)
df <- full_join(param_df, lrt_df, by = group_col) %>%
    select(-c(effect, group, starts_with("term"), npar))

# =====================
# ===== SAVE DATA =====
# =====================

if (opt$unharmonized) {
    filename <- str_glue("sim_test_{opt$feature_roi}_unharmonized.csv")
} else {
    filename <- str_glue("sim_test_{opt$feature_roi}.csv")
}
write_csv(
    df,
    file.path(opt$odir, filename)
)

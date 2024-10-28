# ==============================================================================

# Name: get_sim_results_table.R
# Author: Braden Yang
# Created: 02/14/2024
# Description: create simulation results table (Supplementary Table 2 and 3)

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

source("code/functions.R")

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

mean_sd_str <- function(m, s) {

    # https://stackoverflow.com/questions/34365803/how-to-place-plus-minus-operator-in-text-annotation-of-plot-ggplot2

    s <- if_else(
        is.na(m) & is.na(s),
        NA,
        paste0({{m}}, " \u00B1 ", {{s}})
    )
    return(s)

}

# ============================
# ===== DEFINE VARIABLES =====
# ============================

# harmonization method levels and labels
harm_method_levels <- c(
    "suvr",
    "cl",
    "combat_nocovar",
    "combat_age_sex_apoe",
    "peace__nocovar",
    "peace__age_sex_apoe",
    "longCombat__nocovar",
    "longCombat__age_sex_apoe"
)
harm_method_labels <- c(
    "unharmonized",
    "centiloid",
    "ComBat, no covariates",
    "ComBat + age, sex, APOE",
    "PEACE, no covariates",
    "PEACE\n+ age, sex, APOE",
    "longCombat,\nno covariates",
    "longCombat\n+ age, sex, APOE"
)
harm_method_dict <- harm_method_labels
names(harm_method_dict) <- harm_method_levels

# =====================
# ===== LOAD DATA =====
# =====================

sim_df <- read_rds("data/simulation/stats.RDS")
roi_df <- read_rds("data/csv/roi.RDS")

# =====================
# ===== TIDY DATA =====
# =====================

# summary region
sim_tidy_summary <- sim_df %>%
    filter(
        n_subj == 50,
        treatment.delta %in% c(0, 0.01, 0.02, 0.03),
        feature_roi == "summary"
    ) %>%
    group_by(harmonization_method, treatment.delta) %>%
    summarise(power.m = mean(power), power.sd = sd(power)) %>%
    mutate(
        mean_power_str = mean_sd_str(round(power.m, 3), round(power.sd, 3)),
        treatment.delta = str_c("Annualized rate-of-change = ", -treatment.delta)
    ) %>%
    ungroup() %>%
    select(-starts_with("power")) %>%
    pivot_wider(
        names_from = treatment.delta,
        values_from = mean_power_str
    ) %>%
    harm_method_to_factor() %>%
    arrange(harmonization_method)

# ROIs
sim_tidy_roi <- sim_df %>%
    filter(
        n_subj == 50,
        treatment.delta %in% c(0, 0.01, 0.02, 0.03),
        feature_roi != "summary"
    ) %>%
    left_join(  # add ROI type column
        roi_df %>%
            filter(study == "oasis", feature_type == "suvr") %>%
            select(fs_label, roi_type) %>%
            rename(feature_roi = fs_label),
        by = "feature_roi"
    )
sim_tidy_roi <- sim_tidy_roi %>%
    bind_rows(  # add "All ROI" as an ROI type (compute statistics across all ROI)
        sim_tidy_roi %>% mutate(roi_type = "All ROI")
    )
sim_tidy_roi <- sim_tidy_roi %>%
    group_by(harmonization_method, roi_type, treatment.delta, feature_roi) %>%
    summarise(
        mean_power = mean(power)
    ) %>%
    group_by(harmonization_method, roi_type, treatment.delta) %>%
    summarise(
        mean_power.m = mean(mean_power),
        mean_power.sd = sd(mean_power)
    ) %>%
    mutate(mean_power_str = mean_sd_str(round(mean_power.m, 3), round(mean_power.sd, 3))) %>%
    ungroup()

sim_wide_roi <- sim_tidy_roi %>%
    select(-c(mean_power.m, mean_power.sd)) %>%
    pivot_wider(
        names_from = c(roi_type),
        values_from = mean_power_str
    ) %>%
    harm_method_to_factor() %>%
    arrange(harmonization_method)

# ===========================
# ===== CREATE GT TABLE =====
# ===========================

# summary region
results_summary_gt <- sim_tidy_summary %>%
    gt(rowname_col = "harmonization_method") %>%
    tab_stubhead(label = "Harmonization method") %>%
    tab_spanner(
        label = "Mean power",
        columns = everything()
    )

# ROI
results_roi_gt <- sim_wide_roi %>%
    gt(rowname_col = "harmonization_method") %>%
    tab_stubhead(label = "Harmonization method")

for (t in rev(c(0, -0.01, -0.02, -0.03))) {

    results_roi_gt <- results_roi_gt %>%
        tab_row_group(
            label = str_glue("Annualized rate-of-change = {t}"),
            rows = treatment.delta == -t,
            id = str_glue("treatment.delta.{t}")
        )

}

results_roi_gt <- results_roi_gt %>%
    cols_hide(treatment.delta) %>%
    tab_spanner(
        label = "Mean power",
        columns = everything()
    ) %>%
    cols_move(
        columns = c(`Other cortical ROI`, `Subcortical ROI`, `All ROI`),
        after = `Summary cortical ROI`
    )

# ==============================
# ===== SAVE TABLE TO FILE =====
# ==============================

gtsave(results_summary_gt, "sim_summary_results.html", opt$odir)
gtsave(results_summary_gt, "sim_summary_results.rtf", opt$odir)
gtsave(results_roi_gt, "sim_roi_results.html", opt$odir)
gtsave(results_roi_gt, "sim_roi_results.rtf", opt$odir)

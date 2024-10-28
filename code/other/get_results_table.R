# ==============================================================================

# Name: get_results_table.R
# Author: Braden Yang
# Created: 02/14/2024
# Description: create results table (Table 2)

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

# =====================
# ===== LOAD DATA =====
# =====================

icc_df <- read_rds("data/head2head/icc.RDS")
mae_df <- read_rds("data/head2head/mae.RDS")
icc_ttest <- read_rds("data/head2head/icc_ttest.RDS")
mae_ttest <- read_rds("data/head2head/mae_ttest.RDS")
ae_summary_ttest <- read_rds("data/head2head/ae_summary_ttest.RDS")

# ==========================================
# ===== COMPILE RESULTS INTO DATAFRAME =====
# ==========================================

icc_results <- icc_df %>%
    group_by(harmonization_method) %>%
    summarise(
        icc.summary = icc[feature_roi == "summary"],
        icc.mean_roi = mean(icc[feature_roi != "summary"]),
        icc.mean_roi.sd = sd(icc[feature_roi != "summary"]),
        icc.mean_summary_roi = mean(icc[feature_roi != "summary" & roi_type == "Summary cortical ROI"]),
        icc.mean_summary_roi.sd = sd(icc[feature_roi != "summary" & roi_type == "Summary cortical ROI"]),
        icc.mean_other_roi = mean(icc[feature_roi != "summary" & roi_type == "Other cortical ROI"]),
        icc.mean_other_roi.sd = sd(icc[feature_roi != "summary" & roi_type == "Other cortical ROI"]),
        icc.mean_subctx_roi = mean(icc[feature_roi != "summary" & roi_type == "Subcortical ROI"]),
        icc.mean_subctx_roi.sd = sd(icc[feature_roi != "summary" & roi_type == "Subcortical ROI"])
    ) %>%
    harm_method_to_factor()
icc_allroi_ttest <- icc_ttest %>%
    filter(group1 == "unharmonized") %>%
    select(roi_type, group2, p.adj.signif) %>%
    rename(
        harmonization_method = group2,
        signif = p.adj.signif
    ) %>%
    mutate(
        roi_type = case_match(roi_type,
            "All ROI" ~ "icc.mean_roi",
            "Summary cortical ROI" ~ "icc.mean_summary_roi",
            "Other cortical ROI" ~ "icc.mean_other_roi",
            "Subcortical ROI" ~ "icc.mean_subctx_roi"
        ) %>% str_c(".signif")
    ) %>%
    pivot_wider(
        names_from = roi_type,
        values_from = signif
    )

mae_results <- mae_df %>%
    group_by(harmonization_method) %>%
    summarise(
        mae.summary = mean_abs_diff[feature_roi == "summary"],
        mae.summary.sd = sd_abs_diff[feature_roi == "summary"],
        mae.mean_roi = mean(mean_abs_diff[feature_roi != "summary"]),
        mae.mean_roi.sd = sd(mean_abs_diff[feature_roi != "summary"]),
        mae.mean_summary_roi = mean(mean_abs_diff[feature_roi != "summary" & roi_type == "Summary cortical ROI"]),
        mae.mean_summary_roi.sd = sd(mean_abs_diff[feature_roi != "summary" & roi_type == "Summary cortical ROI"]),
        mae.mean_other_roi = mean(mean_abs_diff[feature_roi != "summary" & roi_type == "Other cortical ROI"]),
        mae.mean_other_roi.sd = sd(mean_abs_diff[feature_roi != "summary" & roi_type == "Other cortical ROI"]),
        mae.mean_subctx_roi = mean(mean_abs_diff[feature_roi != "summary" & roi_type == "Subcortical ROI"]),
        mae.mean_subctx_roi.sd = sd(mean_abs_diff[feature_roi != "summary" & roi_type == "Subcortical ROI"])
    ) %>%
    harm_method_to_factor()
ae_ttest_signif <- ae_summary_ttest %>%
    filter(group1 == "unharmonized") %>%
    select(group2, p.adj.signif) %>%
    rename(
        harmonization_method = group2,
        mae.summary.signif = p.adj.signif
    )
mae_allroi_ttest <- mae_ttest %>%
    filter(group1 == "unharmonized") %>%
    select(roi_type, group2, p.adj.signif) %>%
    rename(
        harmonization_method = group2,
        signif = p.adj.signif
    ) %>%
    mutate(
        roi_type = case_match(roi_type,
            "All ROI" ~ "mae.mean_roi",
            "Summary cortical ROI" ~ "mae.mean_summary_roi",
            "Other cortical ROI" ~ "mae.mean_other_roi",
            "Subcortical ROI" ~ "mae.mean_subctx_roi"
        ) %>% str_c(".signif")
    ) %>%
    pivot_wider(
        names_from = roi_type,
        values_from = signif
    )

results_table <- purrr::reduce(
    list(icc_results, icc_allroi_ttest, mae_results, ae_ttest_signif, mae_allroi_ttest),
    ~full_join(.x, .y, by = "harmonization_method")
) %>%
    mutate(
        across(
            where(is.numeric),
            ~ round(.x, digits = 3)
        ),
        across(
            ends_with(".signif"),
            ~ if_else(.x == "ns" | is.na(.x),"", .x)
        ),
        harmonization_method = factor(harmonization_method, levels = harm_method_labels)
    ) %>%
    arrange(harmonization_method) %>%
    mutate(
        icc.mean_roi = str_glue("{icc.mean_roi} \u00B1 {icc.mean_roi.sd}<sup>{icc.mean_roi.signif}</sup>"),
        icc.mean_summary_roi = str_glue("{icc.mean_summary_roi} \u00B1 {icc.mean_summary_roi.sd}<sup>{icc.mean_summary_roi.signif}</sup>"),
        icc.mean_other_roi = str_glue("{icc.mean_other_roi} \u00B1 {icc.mean_other_roi.sd}<sup>{icc.mean_other_roi.signif}</sup>"),
        icc.mean_subctx_roi = str_glue("{icc.mean_subctx_roi} \u00B1 {icc.mean_subctx_roi.sd}<sup>{icc.mean_subctx_roi.signif}</sup>"),
        mae.summary = str_glue("{mae.summary} \u00B1 {mae.summary.sd}<sup>{mae.summary.signif}</sup>"),
        mae.mean_roi = str_glue("{mae.mean_roi} \u00B1 {mae.mean_roi.sd}<sup>{mae.mean_roi.signif}</sup>"),
        mae.mean_summary_roi = str_glue("{mae.mean_summary_roi} \u00B1 {mae.mean_summary_roi.sd}<sup>{mae.mean_summary_roi.signif}</sup>"),
        mae.mean_other_roi = str_glue("{mae.mean_other_roi} \u00B1 {mae.mean_other_roi.sd}<sup>{mae.mean_other_roi.signif}</sup>"),
        mae.mean_subctx_roi = str_glue("{mae.mean_subctx_roi} \u00B1 {mae.mean_subctx_roi.sd}<sup>{mae.mean_subctx_roi.signif}</sup>"),
        .keep = "unused"
    ) %>%
    mutate(across(starts_with(c("icc", "mae")), ~ map(.x, html)))
# for some reason, the `html` function doesn't work when you apply it to a vector; you need to map it
# individually to each cell in the table, and have a vector of lists which contain html objects. Only then
# does it work uwu

# ===========================
# ===== CREATE GT TABLE =====
# ===========================

results_gt <- results_table %>%
    select(
        harmonization_method,
        icc.summary, mae.summary,
        icc.mean_roi, mae.mean_roi,
        icc.mean_summary_roi, mae.mean_summary_roi,
        icc.mean_other_roi, mae.mean_other_roi,
        icc.mean_subctx_roi, mae.mean_subctx_roi,
    ) %>%
    gt(rowname_col = "harmonization_method") %>%
    tab_stubhead(label = "Harmonization method") %>%
    cols_label(
        icc.summary = "ICC",
        icc.mean_roi = "ICC",
        icc.mean_summary_roi = "ICC",
        icc.mean_other_roi = "ICC",
        icc.mean_subctx_roi = "ICC",
        mae.summary = "MAE",
        mae.mean_roi = "MAE",
        mae.mean_summary_roi = "MAE",
        mae.mean_other_roi = "MAE",
        mae.mean_subctx_roi = "MAE",
    ) %>%
    tab_spanner(
        label = "Global summary",
        columns = c("icc.summary", "mae.summary"),
        level = 2
    ) %>%
    tab_spanner(
        label = "All ROI",
        columns = c("icc.mean_roi", "mae.mean_roi")
    ) %>%
    tab_spanner(
        label = "Summary cortical ROI",
        columns = c("icc.mean_summary_roi", "mae.mean_summary_roi")
    ) %>%
    tab_spanner(
        label = "Other cortical ROI",
        columns = c("icc.mean_other_roi", "mae.mean_other_roi")
    ) %>%
    tab_spanner(
        label = "Subcortical ROI",
        columns = c("icc.mean_subctx_roi", "mae.mean_subctx_roi")
    ) %>%
    tab_spanner(
        label = "FreeSurfer ROI",
        spanners = c("All ROI", "Summary cortical ROI", "Other cortical ROI", "Subcortical ROI")
    )

# ==============================
# ===== SAVE TABLE TO FILE =====
# ==============================

gtsave(results_gt, "results.html", opt$odir)
gtsave(results_gt, "results.rtf", opt$odir)

# ==============================================================================

# Name: crossover_analysis.R
# Author: Braden Yang
# Created: 06/30/2023
# Description: evaluate harmonization of OASIS crossover dataset, including ICC
#   analysis, MAE analysis, and paired t-tests

# ==============================================================================

rm(list = ls())

# *** toggle variables ***
INTERACTIVE <- FALSE
SAVE_FIG <- TRUE
# *** toggle variables ***

# ===========================
# ===== IMPORT PACKAGES =====
# ===========================

library(tidyverse)
library(optparse)

library(rstatix)
library(psych)

library(patchwork)
library(ggpubr)
library(ggseg)

library(devtools)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store",
        type="character", help="Path to project directory (if none, uses current working directory)"),
    make_option(c("-o", "--odir"), action="store", default="figures/crossover",
        type="character", help="Path to output directory"),
    make_option(c("-p", "--pvc"), action="store_true", default=FALSE,
        help="if specified, run analyses using PVC data"),
    make_option(c("-c", "--combat_dir", action="store", default=NULL,
        help="path to alternate ComBat directory")),
    make_option(c("-e", "--plot_peace"), action="store_true", default=FALSE,
        help="plot PEACE results")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (INTERACTIVE) {
    opt$plot_peace <- FALSE
}

# set working directory
if (!is.null(opt$wdir)) {
    setwd(opt$wdir)
}

# load byyfunctions
load_all("submodules/byyfunctions")

# ======================================
# ===== SOURCE FUNCTIONS/VARIABLES =====
# ======================================

source("code/functions.R")

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

icc_single_fixed <- function(.data) {

    # compute ICC of a single fixed model using `psych` R library

    icc_result <- ICC(.data, lmer = FALSE)
    return(icc_result$results[[3, 2]]) # pick out single fixed ICC

}

get_icc_label <- function(icc) {
    return(
        paste0("ICC = ", signif(icc, digits = 3))
    )
}

get_roi_subgroup <- function(.data) {
    .data <- .data %>% mutate(  # categorize regions into subgroups
        roi_type = case_when(
            feature_roi %in% roi_subgroup_list[[1]] ~ "Summary cortical ROI",
            feature_roi %in% roi_subgroup_list[[2]] ~ "Other cortical ROI",
            feature_roi %in% roi_subgroup_list[[3]] ~ "Subcortical ROI",
            .default = NA
        ) %>% factor(levels = c("Summary cortical ROI", "Other cortical ROI", "Subcortical ROI"))
    )
    return(.data)
}

# ============================
# ===== DEFINE VARIABLES =====
# ============================

# get list of ROIs
feature_names <- read_csv("data/csv/feature_names.csv", col_names = FALSE) %>% pull(1)

# get ROI sub-groups
roi_df <- read_rds("data/csv/roi.RDS")
roi_subgroup_names <- c("Summary cortical ROI", "Other cortical ROI", "Subcortical ROI")
roi_subgroup_list <- map(
    roi_subgroup_names,
    ~ roi_df %>%
        filter(study == "oasis", feature_type == "suvr") %>%
        filter(roi_type == .x) %>%
        pull(fs_label),
)

# define output directories
icc_odir <- file.path(opt$odir, "icc")
icc_surf_odir <- file.path(icc_odir, "surf")
mae_odir <- file.path(opt$odir, "mae")
mae_surf_odir <- file.path(mae_odir, "surf")
peace_odir <- file.path(opt$odir, "PEACE")

byyfunctions::make_dir(icc_surf_odir)
byyfunctions::make_dir(mae_surf_odir)
if (opt$plot_peace) {byyfunctions::make_dir(peace_odir)}

# =====================
# ===== LOAD DATA =====
# =====================

if (opt$pvc) {
    suvr_path <- "data/crossover_pvc.RDS"
    cl_path <- "data/centiloid/centiloid_crossover_pvc.RDS"
    if (!is.null(opt$combat_dir)) {
        combat_dir <- opt$combat_dir
    } else {
        combat_dir <- "data/combat/pvc"
    }
} else {
    suvr_path <- "data/crossover.RDS"
    cl_path <- "data/centiloid/centiloid_crossover.RDS"
    if (!is.null(opt$combat_dir)) {
        combat_dir <- opt$combat_dir
    } else {
        combat_dir <- "data/combat"
    }
}

# load crossover data
# NOTE: round age to nearest decimal to address RDS vs. CSV float value differences
crossover_list <- list()
crossover_list[["suvr"]] <- read_rds(suvr_path) %>% mutate(across(starts_with("age"), round_age))
crossover_list[["cl"]] <- read_rds(cl_path) %>% mutate(across(starts_with("age"), round_age))
crossover_df <- bind_rows(crossover_list, .id = "harmonization_method")

# load all combat single covariate experimental results
crossover_combat_files <- list.files(path = combat_dir, pattern = "*.csv")
crossover_combat_df <- read_csv(file.path(combat_dir, crossover_combat_files), id = "file", show_col_types = FALSE) %>%
    mutate(
        covariate = file %>% basename %>% str_remove(".csv") %>% str_remove("combat_crossover_"),
        harmonization_method = paste0(covariate)
    ) %>%
    select(-c(covariate, file)) %>%
    mutate_combat()

# REVISED: filter out crossover data that doesn't appear in ComBat outputs (due to NA values in batch effect variable)
combat_idx_filter <- crossover_combat_df %>% pull(idx) %>% unique

# load PEACE results
if (opt$plot_peace) {
    peace_files <- list.files(path = "data/head2head/PEACE", pattern = "crossover_PEACE.*\\.RDS")
    peace_list <- list()
    for (i in 1:length(peace_files)) {
        peace_list[[i]] <- read_rds(file.path("data/head2head/PEACE", peace_files[[i]])) %>% mutate_combat()
    }
    peace_df <- purrr::reduce(peace_list, bind_rows)
}

# join with main crossover dataframe
if (opt$plot_peace) {
    crossover_df <- bind_rows(
        crossover_df %>% filter(idx %in% combat_idx_filter),
        crossover_combat_df,
        peace_df
    )
} else {
    crossover_df <- bind_rows(
        crossover_df %>% filter(idx %in% combat_idx_filter),
        crossover_combat_df
    )
}

# separate orig methods with new (revisions) methods
orig_harm_methods <- harm_method_levels[!(harm_method_levels %>% str_detect("peace"))]
peace_harm_methods <- harm_method_levels[(harm_method_levels %>% str_detect("peace"))]

# =====================
# ===== TIDY DATA =====
# =====================

crossover_long <- crossover_df %>%
    select(
        subj, pair_id, tracer, harmonization_method,
        sex, apoe,
        ends_with(c(".AV45", ".PIB")),
        all_of(feature_names)) %>% # select only relevant cols (ensures proper pivoting)
    pivot_longer(
        cols = all_of(feature_names),
        names_to = "feature_roi",
        values_to = "feature_value"
    ) %>%
    pivot_wider(
        names_from = "tracer",
        values_from = "feature_value"
    ) %>%
    get_roi_subgroup() %>%
    filter(!(feature_roi %in% c("left.cerebellum.cortex", "right.cerebellum.cortex")))  # exclude cerebellum cortex

# REVISIONS: get corresponding SUVRs at each anchor point (0CL, 100CL)
cl_to_suvr <- function(cl) {
    av45_0cl <- 181 / 163.6
    av45_100cl <- (100 + 181) / 163.6
    pib_0cl <- 119.3 / 111.8
    pib_100cl <- (100 + 119.3) / 111.8
    avg_0cl <- (av45_0cl + pib_0cl) / 2
    avg_100cl <- (av45_100cl + pib_100cl) / 2

    suvr <- cl * (avg_100cl - avg_0cl) / 100 + avg_0cl
    return(suvr)
}
crossover_long <- crossover_long %>%
    mutate(across(
        c(AV45, PIB),
        ~ if_else(
            harmonization_method == "cl",
            cl_to_suvr(.x),
            .x
        )
    ))

# ======================================================
# ===== COMPUTE INTRACLASS CORRELATION COEFFICIENT =====
# ======================================================

# compute ICC
# https://stackoverflow.com/questions/62182231/applying-simple-function-via-across-within-nested-data-on-each-group
# https://stackoverflow.com/questions/75670222/mutate-across-pairs-of-columns-with-dplyr
icc_long <- crossover_long %>%
    nest(.by = c("feature_roi", "roi_type", "harmonization_method")) %>%
    rowwise() %>%
    mutate(icc = data.frame(data$AV45, data$PIB) %>% icc_single_fixed()) %>%
    ungroup()

# compute difference in ICC
icc_diff <- icc_long %>%
    select(-data) %>%
    byyfunctions::pivot_difference(
        names_col = "harmonization_method",
        values_col = "icc",
        ref = "suvr"
    )

# =======================================
# ===== COMPUTE MEAN ABSOLUTE ERROR =====
# =======================================

# compute absolute difference
abs_diff_df <- crossover_long %>%
    mutate(abs_diff = abs(AV45 - PIB))

# summarise mean absolute error
mae_df <- abs_diff_df %>%
    group_by(harmonization_method, feature_roi, roi_type) %>%
    summarise(
        mean_abs_diff = mean(abs_diff),
        sd_abs_diff = sd(abs_diff),
        .groups = "drop"
    )

# compute difference in MAE
mae_diff <- mae_df %>%
    select(-sd_abs_diff) %>%
    byyfunctions::pivot_difference(
        names_col = "harmonization_method",
        values_col = "mean_abs_diff",
        ref = "suvr"
    ) %>%
    left_join(
        mae_df %>% select(feature_roi, harmonization_method, sd_abs_diff),
        by = c("feature_roi", "harmonization_method")
    )

# ===========================================
# ===== FIGURES FOR SUMMARY REGION ONLY =====
# ===========================================

# ICC of summary region
plot_icc_scatter <- function(.data, h) {

    data_df <- .data %>%
        filter(
            feature_roi == "summary",
            harmonization_method == h
        )
    min_lim <- min(c(data_df %>% unnest(data) %>% pull(AV45), data_df %>% unnest(data) %>% pull(PIB))) * 0.99
    max_lim <- max(c(data_df %>% unnest(data) %>% pull(AV45), data_df %>% unnest(data) %>% pull(PIB))) * 1.01
    icc_str <- get_icc_label(data_df$icc)
    p <- data_df %>%
        unnest(data) %>%
        harm_method_to_factor() %>%
        {ggplot(data = ., aes(x = AV45, y = PIB)) +
            geom_point(size = 0.5) +
            geom_smooth(aes(color = harmonization_method), method = "lm") + 
            geom_abline(slope = 1, intercept = 0, color = "#575757") +
            geom_text( # add labels to indicate ICC value
                mapping = aes(label = icc_str, x = Inf, y = -Inf, hjust = 1, vjust = -0.4),
                size = 5
            ) +
            theme_classic() +
            theme(
                legend.position = "none",
                axis.title = element_blank(),
                plot.title = element_text(hjust = 0.5, size = 18),
                axis.text = element_text(size = 12)
            ) +
            coord_fixed(ratio = 1, xlim = c(min_lim, max_lim), ylim = c(min_lim, max_lim)) +
            labs(title = .$harmonization_method %>% unique %>% str_wrap(20))}
    
    return(p)

}

summary_icc_scatter_list <- map(orig_harm_methods, ~ icc_long %>% plot_icc_scatter(.x))
if (opt$plot_peace) {
    peace_summary_icc_scatter_list <- map(peace_harm_methods, ~ icc_long %>% plot_icc_scatter(.x))
    summary_icc_scatter_list <- c(summary_icc_scatter_list, peace_summary_icc_scatter_list)
    x_label_idx <- 4
} else {
    x_label_idx <- 3
}
summary_icc_scatter_list[[1]] <- summary_icc_scatter_list[[1]] + labs(y = "PiB measurement") + theme(axis.title.y = element_text(size = 20))
summary_icc_scatter_list[[x_label_idx]] <- summary_icc_scatter_list[[x_label_idx]] + labs(x = "FBP measurement") + theme(axis.title.x = element_text(size = 20))
summary_icc_scatter <- wrap_plots(summary_icc_scatter_list, nrow = 1)

# absolute difference of summary region across scans
abs_diff_summary_paired_ttest <- TRUE  # whether to perform paired t-test
abs_diff_summary_ttest <- abs_diff_df %>%
    filter(feature_roi == "summary") %>%
    harm_method_to_factor() %>%
    t_test(
        abs_diff ~ harmonization_method,
        p.adjust.method = "bonferroni",
        paired = abs_diff_summary_paired_ttest
    ) %>%
    add_xy_position(
        x = "harmonization_method",
        fun = "max"
    )

if (opt$plot_peace) {
    plot_df <- abs_diff_df
    pval_df <- abs_diff_summary_ttest
} else {
    plot_df <- abs_diff_df %>%
        filter(!str_detect(harmonization_method, "peace"))
    pval_df <- abs_diff_summary_ttest %>%
        filter(!(str_detect(group1, "PEACE") | str_detect(group2, "PEACE")))
}

abs_diff_summary_boxplot <- plot_df %>%
    filter(feature_roi == "summary") %>%
    harm_method_to_factor() %>%
    ggplot(aes(x = harmonization_method, y = abs_diff, color = harmonization_method)) +
        geom_boxplot(fill = NA, outlier.shape = NA, width = 0.75) +
        geom_point(fill = NA, position = position_jitterdodge(jitter.width = 1.5), shape = 1, size = 0.75) +
        labs(x = "harmonization method", y = "absolute error") +
        theme_classic() +
        theme(axis.text.x = element_text(hjust = 0.5, vjust = 1), legend.position = "none") +
        stat_pvalue_manual(
            pval_df,
            label = "p.adj.signif",
            hide.ns = TRUE,
            tip.length = 0.01,
            step.increase = 0.025,
            bracket.nudge.y = 0.01
        )

abs_diff_summary_boxplot_horiz <- plot_df %>%
    filter(feature_roi == "summary") %>%
    harm_method_to_factor() %>%
    ggplot(aes(x = harmonization_method, y = abs_diff, color = harmonization_method)) +
        geom_boxplot(fill = NA, outlier.shape = NA, width = 0.75) +
        geom_point(fill = NA, position = position_jitterdodge(jitter.width = 1.5), shape = 1, size = 0.75) +
        labs(x = "harmonization method", y = "absolute error") +
        theme_classic() +
        theme(
            axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 16),
            legend.position = "none"
        ) +
        stat_pvalue_manual(
            pval_df,
            label = "p.adj.signif",
            hide.ns = TRUE,
            tip.length = 0.01,
            step.increase = 0.025,
            bracket.nudge.y = 0.01,
            coord.flip = TRUE
        ) +
        coord_flip()

# save figures
if (SAVE_FIG && !opt$plot_peace) {
    ggsave(
        file.path(icc_odir, "scatter_summary.png"),
        summary_icc_scatter,
        width = 16, height = 4, units = "in",
        dpi = 500
    )
    ggsave(
        file.path(mae_odir, "summary_boxplot_horiz.png"),
        abs_diff_summary_boxplot_horiz,
        width = 10, height = 4, units = "in",
        dpi = 500
    )
}

# save PEACE figures (only run if opt$plot_peace==TRUE)
if (SAVE_FIG && opt$plot_peace) {
    ggsave(
        file.path(peace_odir, "scatter_summary_PEACE.png"),
        summary_icc_scatter,
        width = 21, height = 4, units = "in",
        dpi = 500
    )
    ggsave(
        file.path(peace_odir, "summary_boxplot_horiz_PEACE.png"),
        abs_diff_summary_boxplot_horiz,
        width = 15, height = 6, units = "in",
        dpi = 500
    )
}

# ====================================================
# ===== BOXPLOT OF ICC DISTRIUBTIONS ACROSS ROIS =====
# ====================================================

# include roi_type of "All ROI"
icc_long_all <- icc_long %>%
    mutate(roi_type = "All ROI") %>%
    bind_rows(icc_long) %>%
    mutate(roi_type = fct(roi_type, levels = c("All ROI", "Summary cortical ROI", "Other cortical ROI", "Subcortical ROI")))

# perform 2-sample t-test
# https://www.datanovia.com/en/blog/how-to-add-p-values-onto-a-grouped-ggplot-using-the-ggpubr-r-package/
icc_paired <- TRUE  # whether to perform paired t-test
icc_ttest <- icc_long_all %>%
    filter(
        feature_roi != "summary",
    ) %>%
    harm_method_to_factor() %>%
    group_by(roi_type) %>%
    t_test(
        icc ~ harmonization_method,
        p.adjust.method = "bonferroni",
        paired = icc_paired
    ) %>%
    ungroup() %>%
    mutate(
        p.adj = p.adjust(p, method = "bonferroni")
    ) %>%
    add_significance("p.adj")

# plot horizontally and with facets
icc_ttest_horiz <- icc_ttest %>%
    add_xy_position(
        x = "harmonization_method",
        fun = "max"
    ) %>%
    mutate(  # here we manually offset y-position for each facet, since idk how to do it within the ggpubr functions
        y.position = case_match(roi_type,
            "Summary cortical ROI" ~ y.position - 0.04,
            "Other cortical ROI" ~ y.position - 0.1,
            .default = y.position
        )
    )
icc_grouped_boxplot_horiz <- icc_long_all %>%
    filter(feature_roi != "summary") %>%
    harm_method_to_factor() %>%
    ggplot(aes(x = harmonization_method, y = icc, color = harmonization_method)) +
        geom_boxplot(fill = NA, outlier.shape = NA, width = 0.75) +
        geom_point(fill = NA, position = position_jitterdodge(jitter.width = 1.5), shape = 1, size = 0.5) +
        # geom_boxplot(fill = NA, outlier.size = 1, outlier.fill = NA, width = 0.75) +
        facet_wrap(vars(roi_type), ncol = 2, scales = "free") +
        guides(color = "none") +
        labs(x = "ROI subgroup", y = "intraclass correlation coefficient (ICC)") +
        theme_classic() + 
        theme(
            axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 12),
            # legend.position.inside = c(0.85, 0.85),
            panel.spacing.y = unit(0.25, "in"),
            panel.spacing.x = unit(0.4, "in"),
            strip.background = element_blank(),
            strip.text = element_text(size = 14, face = "bold", family = "Arial"),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_blank()
        ) + 
        stat_pvalue_manual(
            icc_ttest_horiz,
            label = "p.adj.signif",
            hide.ns = TRUE,
            tip.length = 0.01,
            step.increase = 0.025,
            bracket.nudge.y = 0.01,
            coord.flip = TRUE
        ) +
        coord_flip()

if (SAVE_FIG && !opt$plot_peace) {
    ggsave(
        file.path(icc_odir, "boxplot_horiz.png"),
        icc_grouped_boxplot_horiz,
        width = 14, height = 5, units = "in",
        dpi = 500
    )
}
if (SAVE_FIG && opt$plot_peace) {
    ggsave(
        file.path(peace_odir, "icc_boxplot_horiz_PEACE.png"),
        icc_grouped_boxplot_horiz,
        width = 18, height = 6, units = "in",
        dpi = 500
    )
}

# ==========================
# ===== ICC ON SURFACE =====
# ==========================

# ICC on surface
# plot for each harmonization method individually
for (h in orig_harm_methods) {

    icc_surf_scale <- scale_fill_viridis_c(
        limits = c(
            icc_long %>% filter(harmonization_method == h) %>% pull(icc) %>% min,
            icc_long %>% filter(harmonization_method == h) %>% pull(icc) %>% max
        )
    )

    icc_surf_aseg <- icc_long %>%
        filter(harmonization_method == h) %>%
        harm_method_to_factor() %>%
        byyfunctions::add_ggseg_label("feature_roi") %>%
        byyfunctions::plot_surf_aseg("icc", harm_method_dict[[h]], cb_title = "ICC", aseg_side = "coronal") &
            icc_surf_scale
        
    if (SAVE_FIG) {
        ggsave(
            file.path(opt$odir, "icc", "surf", paste0(h, ".png")),
            icc_surf_aseg,
            width = 4, height = 10, units = "in"
        )
    }

}

# ICC difference of harmonization methods
if (opt$plot_peace) {
    hm_filter <- icc_diff %>% pull(harmonization_method) %>% unique
} else {
    hm_filter <- orig_harm_methods
}
icc_diff_scale <- scale_fill_gradient2(
    low = "blue4",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(
        min(icc_diff %>% filter(harmonization_method %in% hm_filter) %>% pull("icc.diff"), na.rm = TRUE),
        max(icc_diff %>% filter(harmonization_method %in% hm_filter) %>% pull("icc.diff"), na.rm = TRUE)
    )
)
icc_diff_harm_surf_list <- icc_diff %>%
    filter(harmonization_method != "suvr", harmonization_method %in% hm_filter) %>%
    byyfunctions::add_ggseg_label("feature_roi") %>%
    harm_method_to_factor() %>%
    byyfunctions::plot_ggseg_patchwork(
        fill_col = "icc.diff",
        group_col = "harmonization_method",
        scale = icc_diff_scale,
        cb_title = "difference in ICC",
        aseg_side = "coronal"
    )

if (opt$plot_peace) {
    ggseg_layout_peace <- "
ABCDEF
GHIJKL
"
    icc_diff_harm_surf <- icc_diff_harm_surf_list %>%
        wrap_plots(design = ggseg_layout_peace, guides = "collect")
} else {
    icc_diff_harm_surf <- icc_diff_harm_surf_list %>%
        wrap_plots(design = ggseg_layout, guides = "collect")
}

# save figures
if (SAVE_FIG && !opt$plot_peace) {
    ggsave(
        file.path(icc_surf_odir, "icc_diff_harm_surf.png"),
        icc_diff_harm_surf,
        width = 10, height = 4, units = "in"
    )
}

if (SAVE_FIG && opt$plot_peace) {
    ggsave(
        file.path(peace_odir, "icc_diff_harm_surf_PEACE.png"),
        icc_diff_harm_surf,
        width = 15, height = 4, units = "in"
    )
}

# ====================================================
# ===== BOXPLOT OF MAE DISTRIUBTIONS ACROSS ROIS =====
# ====================================================

# include roi_type of "All ROI"
mae_df_all <- mae_df %>%
    mutate(roi_type = "All ROI") %>%
    bind_rows(mae_df) %>%
    mutate(roi_type = fct(roi_type, levels = c("All ROI", "Summary cortical ROI", "Other cortical ROI", "Subcortical ROI")))

# perform 2-sample t-test
# https://www.datanovia.com/en/blog/how-to-add-p-values-onto-a-grouped-ggplot-using-the-ggpubr-r-package/
mae_paired <- TRUE  # whether to perform paired t-test
mae_ttest <- mae_df_all %>%
    filter(feature_roi != "summary",) %>%
    harm_method_to_factor() %>%
    group_by(roi_type) %>%
    t_test(
        mean_abs_diff ~ harmonization_method,
        p.adjust.method = "bonferroni",
        paired = mae_paired
    ) %>%
    ungroup() %>%
    mutate(
        p.adj = p.adjust(p, method = "bonferroni")
    ) %>%
    add_significance("p.adj")

# plot horizontally and with facets
mae_ttest_horiz <- mae_ttest %>%
    add_xy_position(
        x = "harmonization_method",
        fun = "max"
    )

if (opt$plot_peace) {
    plot_df <- mae_df_all
    pval_df <- mae_ttest_horiz
} else {
    plot_df <- mae_df_all %>%
        filter(!str_detect(harmonization_method, "peace"))
    pval_df <- mae_ttest_horiz %>%
        filter(!(str_detect(group1, "PEACE") | str_detect(group2, "PEACE")))
}

mae_grouped_boxplot_horiz <- plot_df %>%
    filter(feature_roi != "summary") %>%
    harm_method_to_factor() %>%
    ggplot(aes(x = harmonization_method, y = mean_abs_diff, color = harmonization_method)) +
        geom_boxplot(fill = NA, outlier.shape = NA, width = 0.75) +
        geom_point(fill = NA, position = position_jitterdodge(jitter.width = 1.5), shape = 1, size = 0.5) +
        # geom_boxplot(fill = NA, outlier.size = 1, outlier.fill = NA, width = 0.75) +
        facet_wrap(vars(roi_type), ncol = if (opt$plot_peace) {2} else {2}, scales = "free") +
        guides(color = "none") +
        labs(x = "ROI subgroup", y = "mean absolute error (MAE)") +
        theme_classic() + 
        theme(
            axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 12),
            # legend.position.inside = c(0.85, 0.85),
            panel.spacing.y = unit(0.25, "in"),
            panel.spacing.x = unit(0.4, "in"),
            strip.background = element_blank(),
            strip.text = element_text(size = 14, face = "bold", family = "Arial"),
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_blank(),
            plot.margin = margin(r = 10, unit = "pt")
        ) + 
        stat_pvalue_manual(
            pval_df,
            label = "p.adj.signif",
            hide.ns = TRUE,
            tip.length = 0.01,
            coord.flip = TRUE
        ) +
        coord_flip()

if (SAVE_FIG && !opt$plot_peace) {
    ggsave(
        file.path(mae_odir, "boxplot_horiz.png"),
        mae_grouped_boxplot_horiz,
        width = 14, height = 5, units = "in",
        dpi = 500
    )
}
if (SAVE_FIG && opt$plot_peace) {
    ggsave(
        file.path(peace_odir, "mae_roi_boxplot_horiz_PEACE.png"),
        mae_grouped_boxplot_horiz,
        width = 18, height = 6, units = "in",
        dpi = 500
    )
}

# ==========================
# ===== MAE ON SURFACE =====
# ==========================

# MAE on surface
# plot for each harmonization method individually
for (h in orig_harm_methods) {

    mae_scale <- scale_fill_viridis_c(
        limits = c(
            mae_df %>% filter(harmonization_method == h) %>% pull(mean_abs_diff) %>% min,
            mae_df %>% filter(harmonization_method == h) %>% pull(mean_abs_diff) %>% max
        )
    )

    mae_surf_aseg <- mae_df %>%
        filter(harmonization_method == h) %>%
        harm_method_to_factor() %>%
        byyfunctions::add_ggseg_label("feature_roi") %>%
        byyfunctions::plot_surf_aseg("mean_abs_diff", harm_method_dict[[h]], cb_title = "MAE", aseg_side = "coronal") &
            mae_scale
        
    if (SAVE_FIG) {
        ggsave(
            file.path(opt$odir, "mae", "surf", paste0(h, ".png")),
            mae_surf_aseg,
            width = 4, height = 10, units = "in"
        )
    }

}


# ICC difference of harmonization methods
if (opt$plot_peace) {
    hm_filter <- icc_diff %>% pull(harmonization_method) %>% unique
} else {
    hm_filter <- orig_harm_methods
}
mae_diff_scale <- scale_fill_gradient2(
    low = "blue4",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(
        min(mae_diff %>% filter(harmonization_method %in% hm_filter) %>% pull("mean_abs_diff.diff"), na.rm = TRUE), 
        max(mae_diff %>% filter(harmonization_method %in% hm_filter) %>% pull("mean_abs_diff.diff"), na.rm = TRUE)
    )
)
mae_diff_harm_surf_list <- mae_diff %>%
    filter(harmonization_method != "suvr", harmonization_method %in% hm_filter) %>%
    byyfunctions::add_ggseg_label("feature_roi") %>%
    harm_method_to_factor() %>%
    byyfunctions::plot_ggseg_patchwork(
        fill_col = "mean_abs_diff.diff",
        group_col = "harmonization_method",
        scale = mae_diff_scale,
        cb_title = "difference in MAE",
        aseg_side = "coronal"
    )

if (opt$plot_peace) {
    ggseg_layout_peace <- "
ABCDEF
GHIJKL
"
    mae_diff_harm_surf <- mae_diff_harm_surf_list %>%
        wrap_plots(design = ggseg_layout_peace, guides = "collect")
} else {
    mae_diff_harm_surf <- mae_diff_harm_surf_list %>%
        wrap_plots(design = ggseg_layout, guides = "collect")
}

# save figures
if (SAVE_FIG && !opt$plot_peace) {
    ggsave(
        file.path(mae_surf_odir, "mae_diff_harm_surf.png"),
        mae_diff_harm_surf,
        width = 10, height = 4, units = "in"
    )
}
if (SAVE_FIG && opt$plot_peace) {
    ggsave(
        file.path(peace_odir, "mae_diff_harm_surf_PEACE.png"),
        mae_diff_harm_surf,
        width = 15, height = 4, units = "in"
    )
}

# =======================
# ===== SAVE TABLES =====
# =======================

# odir_head2head <- file.path("data/head2head")
# byyfunctions::make_dir(odir_head2head)

# write_rds(icc_diff, file.path(odir_head2head, "icc.RDS"))
# write_rds(mae_diff, file.path(odir_head2head, "mae.RDS"))
# write_rds(icc_ttest, file.path(odir_head2head, "icc_ttest.RDS"))
# write_rds(mae_ttest, file.path(odir_head2head, "mae_ttest.RDS"))
# write_rds(abs_diff_summary_ttest, file.path(odir_head2head, "ae_summary_ttest.RDS"))

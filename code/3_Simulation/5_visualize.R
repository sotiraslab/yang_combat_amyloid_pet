# ==============================================================================

# Name: 5_visualize.R
# Author: Braden Yang
# Created: 12/26/2023
# Description: Load all simulation results and visualize

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
library(patchwork)
library(ggseg)

library(devtools)

# ===========================
# ===== PARSE ARGUMENTS =====
# ===========================

option_list <- list(
    make_option(c("-w", "--wdir"), action="store", default=NULL,
        type="character", help="Path to project directory (if none, uses current working directory)"),
    make_option(c("-s", "--stats"), action="store", default=NULL,
        type="character", help="path to `stats.RDS`; generated when --save_rds is specified"),
    make_option(c("-o", "--odir"), action="store", default=FALSE, type="character",
        help="path to output directory"),
    make_option(c("-r", "--save_rds"), action="store_true", default=FALSE,
        help="if specified, load all simulation results, store in a single table and save as CSV"),
    make_option(c("-d", "--sim_dir"), action="store", default=NULL, type="character",
        help="path to directory containing simulation outputs; only used if --save_rds is specified"),
    make_option(c("-t", "--other_models"), action="store_true", default=FALSE,
        help="include PEACE and longitudinal ComBat in figures")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (INTERACTIVE) {

    opt$save_rds <- FALSE
    opt$stat_rds <- "data/simulation/stats.RDS"
    
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

load_stats_all <- function(d) {

    load_stats <- function(f, alpha = 0.05) {

        stats_df <-
            read_csv(f, show_col_types = FALSE) %>%
            group_by(harmonization_method, feature_roi) %>%
            mutate(
                # p.sig = `Pr.Chi.` < alpha
                p.sig = sig_one_sided(`Pr.Chi.`, estimate, alpha = alpha, direction_negative = TRUE)  # compute one-sided p-value
            ) %>%
            summarise(
                power = sum(p.sig) / n(),
                estimate.m = mean(estimate),
                estimate.sd = sd(estimate),
                .groups = "drop"
            )
        
        return(stats_df)

    }

    csv_list <- list.files(d, "sim_test_*", full.names = TRUE)
    stats_df <- map(csv_list, load_stats) %>% bind_rows()

    return(stats_df)

}

param_to_factor <- function(.data) {

    n_subj_v <- .data$n_subj %>% unique %>% sort
    p_av45_v <- .data$p_av45.placebo %>% unique %>% sort
    treatment.delta_v <- .data$treatment.delta %>% unique %>% sort

    df <- .data %>%
        mutate(
            n_subj.fct = factor(n_subj, levels = n_subj_v, labels = str_c(n_subj_v, " subjects per group")),
            across(
                starts_with("p_av45"),
                ~ factor(.x, levels = p_av45_v, labels = str_c(floor(p_av45_v*100), "%")),
                .names = "{.col}.fct"
            ),
            treatment.delta.fct = factor(treatment.delta, levels = treatment.delta_v, labels = str_wrap(str_c("annualized rate-of-change = ", treatment.delta_v * -1), width = 25)),
            harmonization_method.fct = factor(harmonization_method, levels = harm_method_levels, labels = harm_method_labels)
        )
    
    return(df)

}

plot_heatmap <- function(
    .data,
    metric_name,
    cb_title,
    cb_height = 5,
    wrap_length = 20,
    row_facet = "harmonization_method.fct",
    col_facet = "treatment.delta.fct",
    no_labels = FALSE,
    no_xticks = FALSE,
    plot_margin = NULL
) {

    # wrapped text: https://stackoverflow.com/questions/21878974/wrap-long-axis-labels-via-labeller-label-wrap-in-ggplot2

    p <- .data %>% ggplot(aes(x = p_av45.placebo.fct, y = p_av45.treatment.fct, fill = .data[[metric_name]])) +
        geom_tile(color = "gray75") +
        facet_grid(rows = vars(.data[[row_facet]]), cols = vars(.data[[col_facet]])) +
        coord_fixed() +
        labs(
            x = if (no_labels) {NULL} else {"% FBP of placebo group"},
            y = if (no_labels) {NULL} else {"% FBP of treatment group"}
        ) +
        guides(fill = guide_colorbar(title = str_wrap(cb_title, width = wrap_length), barheight = cb_height)) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()
        )

        if (!is.null(plot_margin)) {
            p <- p + theme(plot.margin = plot_margin)
        }

        if (no_xticks) {
            p <- p + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
        }

    return(p)

}

plot_diff_heatmap <- function(
    .data,
    metric_name,
    diff_name,
    prefix,
    title = NULL,
    row_facet = "harmonization_method.fct",
    h = NULL,
    return_separate = FALSE
) {

    # plot unharmonized
    sim_heatmap_unharm <- .data %>%
        filter(harmonization_method == "suvr") %>%
        {plot_heatmap(
            .,
            metric_name,
            prefix,
            row_facet = row_facet,
            cb_height = 7,
            no_labels = TRUE,
            no_xticks = TRUE,
            plot_margin = margin(b = -5)
        ) +
            scale_fill_viridis_c(
                limits = c(
                    min(.[[metric_name]], na.rm = TRUE), 
                    max(.[[metric_name]], na.rm = TRUE)
                )
            )}

    # plot harmonized
    sim_heatmap_harm <- .data %>%
        filter(harmonization_method != "suvr") %>%
        {plot_heatmap(., diff_name, str_c("difference in ", prefix), row_facet = row_facet, cb_height = cb_height_harm, wrap_length = 10) +
            scale_fill_gradient2(
                low = "dodgerblue4",
                mid = "white",
                high = "orangered3",
                midpoint = 0,
                limits = c(
                    min(.[[diff_name]], na.rm = TRUE), 
                    max(.[[diff_name]], na.rm = TRUE)
                )
            ) +
            theme(
                strip.background.x = element_blank(),
                strip.text.x = element_blank(),
                axis.title.y = element_text(hjust = 0.75)
            )}

    # combine
    sim_heatmap <- sim_heatmap_unharm + sim_heatmap_harm +
        plot_layout(heights = h) +
        plot_annotation(
            title = title,
            theme = theme(plot.title = element_text(hjust = 0.5))
        )

    if (return_separate) {
        return(
            list("unharmonized" = sim_heatmap_unharm, "harmonized" = sim_heatmap_harm)
        )
    } else {
        return(sim_heatmap)
    }

}

sig_one_sided <- function(p, effect, alpha = 0.05, direction_negative = FALSE) {

    # convert a 2-sided p-value into a one-sided p-value with a direction
    # - direction_negative = FALSE: positive effect sizes are considered significant
    # - direction_negative = TRUE: negative effect sizes are considered significant
    # https://stats.stackexchange.com/questions/325354/if-and-how-to-use-one-tailed-testing-in-multiple-regression

    if (direction_negative) {
        effect <- effect * -1
    }

    p_half <- p / 2
    p_sig <- if_else(effect >= 0, p_half, 1 - p_half)
    sig <- p_sig < alpha

    return(sig)

}

plot_sim_ggseg <- function(
    .data,
    metric,
    my_scale,
    aseg = FALSE,
    plot_margin = NULL
) {

    if (aseg) {
        g <- ggseg::geom_brain(
            atlas = ggseg::aseg,
            mapping = ggplot2::aes(fill = .data[[metric]]),
            side = "coronal"
        )
    } else {
        g <- ggseg::geom_brain(
            atlas = ggseg::dk,
            position = ggseg::position_brain(hemi ~ side),
            mapping = ggplot2::aes(fill = .data[[metric]])
        )
    }

    p <- .data %>%
        byyfunctions::add_ggseg_label() %>%
        group_by(harmonization_method.fct, treatment.delta.fct) %>%
        ggplot() + g +
            facet_grid(rows = vars(harmonization_method.fct), cols = vars(treatment.delta.fct)) +
            theme(
                panel.background = element_blank(),
                panel.border = element_rect(fill = NA, color = "gray75"),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                plot.title = element_text(hjust = 0.5)
            ) +
            my_scale
    
    if (!is.null(plot_margin)) {
        p <- p + theme(plot.margin = plot_margin)
    }

    return(p)

}

# ============================
# ===== DEFINE VARIABLES =====
# ============================

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
) %>% str_wrap(width = 18)
harm_method_dict <- harm_method_labels
names(harm_method_dict) <- harm_method_levels

if (opt$other_models) {
    cb_height_harm <- 41
} else {
    cb_height_harm <- 23
}

# ==================================================
# ===== LOAD ALL RESULTS AND SAVE AS ONE TABLE =====
# ==================================================

if (opt$save_rds) {

    sim_list <- list.files(opt$sim_dir, "sim_*")
    byyfunctions::make_dir(opt$odir)

    # load all results, compute power/mean effsize, store in single tibble
    sim_stats_df <- tibble(
        d = file.path(opt$sim_dir, sim_list),
        param = sim_list
    )

    sim_stats_df <- sim_stats_df %>%
        tidyr::separate_wider_delim(
            param,
            delim = "_",
            names = c(".rm", "n_subj", "p_av45.placebo", "p_av45.treatment", "treatment.delta")
        ) %>%
        mutate(
            across(
                c(n_subj, starts_with("p_av45"), treatment.delta),
                as.numeric
            )
        ) %>%
        select(-.rm) %>%
        mutate(
            stats = map(d, load_stats_all)
        ) %>%
        unnest(stats)

    # save as RDS
    sim_stats_df %>%
        write_rds(
            file.path(opt$odir, "stats.RDS")
        )

    stop("++ successfully saved output; exiting ++")

}

# ======================================
# ===== OUTPUT FOLDERS FOR FIGURES =====
# ======================================

if (SAVE_FIG) {

    byyfunctions::make_dir(opt$odir)

    odir_surf <- file.path(opt$odir, "surf")
    byyfunctions::make_dir(odir_surf)

}

# ======================
# ===== LOAD STATS =====
# ======================

sim_stats <- read_rds(opt$stats)

if (!opt$other_models) {
    sim_stats <- sim_stats %>% filter(
        harmonization_method %in% c("suvr", "cl", "combat_nocovar", "combat_age_sex_apoe")
    )
}

# compute difference compared to unharmonized
sim_stats_diff <- sim_stats %>%
    byyfunctions::pivot_difference(
        names_col = "harmonization_method",
        values_col = c("power", "estimate.m", "estimate.sd"),
        ref = "suvr"
    )

# filter out effect sizes
sim_stats_diff <- sim_stats_diff %>%
    filter(
        treatment.delta %in% c(0, 0.01, 0.02, 0.03)
    ) %>%
    param_to_factor()

sim_stats_regional_mean_power <- sim_stats_diff %>%
    mutate(rate_of_change.absdiff = abs(- treatment.delta - estimate.m)) %>%
    group_by(pick(starts_with(c("harmonization_method", "n_subj", "treatment.delta"))), feature_roi) %>%
    summarise(power.mean = mean(power), .groups = "drop") %>%
    select(-harmonization_method.fct) %>%
    byyfunctions::pivot_difference(
        names_col = "harmonization_method",
        values_col = "power.mean",
        ref = "suvr"
    ) %>%
    mutate(
        harmonization_method.fct = factor(harmonization_method, levels = harm_method_levels, labels = harm_method_labels)
    )

# =====================================================
# ===== VISUALIZE EFFECT SIZE EXPERIMENT, SUMMARY =====
# =====================================================

power_heatmap <- list()
power_diff_heatmap <- list()
estimate_heatmap <- list()
estimate_diff_heatmap <- list()

roi_filter <- "summary"

for (ns in unique(sim_stats$n_subj)) {

    # filter summary region, n_subj
    sim_stats_filter <- sim_stats_diff %>%
        filter(
            n_subj == ns,
            feature_roi == roi_filter,
            !is.na(harmonization_method.fct)
        )

    # plot change in power relative to SUVR
    # - desired output: increase in power after harmonization, especially in smaller effect sizes
    power_diff_heatmap[[as.character(ns)]] <- sim_stats_filter %>%
        plot_diff_heatmap(
            metric_name = "power",
            diff_name = "power.diff",
            prefix = "power",
            h = c(1, 3),
            return_separate = opt$other_models
        )

    if (SAVE_FIG) {

        if (opt$other_models) {
            ggsave(
                file.path(opt$odir, str_glue("power_diff_{ns}subj_unharmonized.png")),
                power_diff_heatmap[[as.character(ns)]][["unharmonized"]],
                width = 9, height = 3, units = "in", dpi = 300
            )
            ggsave(
                file.path(opt$odir, str_glue("power_diff_{ns}subj_harmonized.png")),
                power_diff_heatmap[[as.character(ns)]][["harmonized"]],
                width = 9, height = 15, units = "in", dpi = 300
            )
        } else {
            ggsave(
                file.path(opt$odir, str_glue("power_diff_{ns}subj.png")),
                power_diff_heatmap[[as.character(ns)]],
                width = 9, height = 8.1, units = "in", dpi = 300
            )
        }

    }

}

# ======================================================
# ===== VISUALIZE EFFECT SIZE EXPERIMENT, REGIONAL =====
# ======================================================

# cortical surface
sim_surf_list <- list()
for (ns in unique(sim_stats$n_subj)) {

    sim_surf_unharm <- sim_stats_regional_mean_power %>%
        filter(harmonization_method == "suvr", n_subj == ns, !is.na(harmonization_method.fct)) %>%
        plot_sim_ggseg("power.mean", scale_fill_viridis_c(), plot_margin = margin(b = -10)) +
            guides(fill = guide_colorbar(title = str_wrap("mean power", width = 15), barheight = 7))
    sim_surf_harm <- sim_stats_regional_mean_power %>%
        filter(harmonization_method != "suvr", n_subj == ns, !is.na(harmonization_method.fct)) %>%
        plot_sim_ggseg(
            "power.mean.diff",
            scale_fill_gradient2(
                low = "dodgerblue4",
                mid = "white",
                high = "orangered3",
                midpoint = 0
            ),
            plot_margin = margin(t = -10)
        ) +
            guides(fill = guide_colorbar(title = str_wrap("difference in mean power", width = 15), barheight = cb_height_harm)) +
            theme(
                strip.background.x = element_blank(),
                strip.text.x = element_blank(),
            )
    
    # combine
    sim_surf_list[[as.character(ns)]] <- sim_surf_unharm + sim_surf_harm +
        plot_layout(heights = c(1,3)) +
        plot_annotation(
            theme = theme(plot.title = element_text(hjust = 0.5))
        )

    if (SAVE_FIG) {
        if (opt$other_models) {
            ggsave(
                file.path(odir_surf, str_glue("power_mean_{ns}subj_unharmonized.png")),
                sim_surf_unharm,
                width = 9, height = 2.5, units = "in", dpi = 300
            )
            ggsave(
                file.path(odir_surf, str_glue("power_mean_{ns}subj_harmonized.png")),
                sim_surf_harm,
                width = 9, height = 10, units = "in", dpi = 300
            )
        } else {
            ggsave(
                file.path(odir_surf, str_glue("power_mean_{ns}subj.png")),
                sim_surf_list[[as.character(ns)]],
                width = 11, height = 8, units = "in", dpi = 300
            )
        }
        
    }
    
}

# subcortical regions
sim_aseg_list <- list()
for (ns in unique(sim_stats$n_subj)) {

    sim_surf_unharm <- sim_stats_regional_mean_power %>%
        filter(harmonization_method == "suvr", n_subj == ns, !is.na(harmonization_method.fct)) %>%
        plot_sim_ggseg("power.mean", scale_fill_viridis_c(), plot_margin = margin(b = -10), aseg = TRUE) +
            guides(fill = guide_colorbar(title = str_wrap("mean power", width = 15), barheight = 7))
    sim_surf_harm <- sim_stats_regional_mean_power %>%
        filter(harmonization_method != "suvr", n_subj == ns, !is.na(harmonization_method.fct)) %>%
        plot_sim_ggseg(
            "power.mean.diff",
            scale_fill_gradient2(
                low = "dodgerblue4",
                mid = "white",
                high = "orangered3",
                midpoint = 0
            ),
            plot_margin = margin(t = -10),
            aseg = TRUE
        ) +
            guides(fill = guide_colorbar(title = str_wrap("difference in mean power", width = 15), barheight = cb_height_harm)) +
            theme(
                strip.background.x = element_blank(),
                strip.text.x = element_blank(),
            )
    
    # combine
    sim_aseg_list[[as.character(ns)]] <- sim_surf_unharm + sim_surf_harm +
        plot_layout(heights = c(1.25,3)) +
        plot_annotation(
            theme = theme(plot.title = element_text(hjust = 0.5))
        )

    if (SAVE_FIG) {
        if (opt$other_models) {
            ggsave(
                file.path(odir_surf, str_glue("aseg_power_mean_{ns}subj_unharmonized.png")),
                sim_surf_unharm,
                width = 9, height = 2.5, units = "in", dpi = 300
            )
            ggsave(
                file.path(odir_surf, str_glue("aseg_power_mean_{ns}subj_harmonized.png")),
                sim_surf_harm,
                width = 9, height = 13, units = "in", dpi = 300
            )
        } else {
            ggsave(
                file.path(odir_surf, str_glue("aseg_power_mean_{ns}subj.png")),
                sim_aseg_list[[as.character(ns)]],
                width = 11, height = 10, units = "in", dpi = 30
            )
        }
    }
    
}

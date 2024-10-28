# ==============================================================================

# Name: get_demographics.R
# Author: Braden Yang
# Created: 02/07/2024
# Description: create demographics table

# ==============================================================================

# clear environment
rm(list = ls())

# ===========================
# ===== IMPORT PACKAGES =====
# ===========================

library(tidyverse)
library(gt)
library(optparse)
library(rstatix)

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

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

select_col <- function(.data) {

    df <- .data %>%
        select(subj, tracer, age, sex, apoe, cdr, amyloid_positive)

    return(df)

}

mean_sd_str <- function(m, s) {

    # https://stackoverflow.com/questions/34365803/how-to-place-plus-minus-operator-in-text-annotation-of-plot-ggplot2

    s <- if_else(
        is.na(m) & is.na(s),
        NA,
        paste0(round({{m}}, 1), " \u00B1 ", round({{s}}, 2))
    )
    return(s)

}

# =======================
# ===== LOAD TABLES =====
# =======================

demo_data_dir <- "data/demographics"
sim_dir <- "data/simulation"

crossover_df <- read_csv(file.path(demo_data_dir, "crossover.csv"))
crossover_train_df <- read_csv(file.path(demo_data_dir, "crossover_train.csv"))
amypos_df <- read_csv(file.path(sim_dir, "amypos.csv"))

# =======================
# ===== TIDY TABLES =====
# =======================

# select: subject, tracer, age at scan, sex, APOE, CDR
data_df <- bind_rows(
    select_col(crossover_df) %>% mutate(cohort = "crossover"),
    select_col(crossover_train_df) %>% mutate(cohort = "crossover_train"),
    select_col(amypos_df) %>% mutate(cohort = "simulation")
)

# compute time interval of each scan
data_df <- data_df %>%
    nest(.by = c("cohort", "tracer")) %>%
    mutate(
        data = map(data, ~ .x %>% byyfunctions::get_diff_adjacent("age", "subj"))
    ) %>%
    unnest(data)

# separate H2H-training and simulation cohorts into single-tracer and mixed
mixed_tracer <- data_df %>%
    select(subj, cohort, tracer) %>%
    filter(cohort %in% c("crossover_train", "simulation")) %>%
    group_by(subj, cohort) %>%
    summarise(subj.mixed.tracer = any(tracer == "AV45") & any(tracer == "PIB"), .groups="drop")
data_df <- data_df %>%
    left_join(mixed_tracer, by = c("subj", "cohort")) %>%
    mutate(
        cohort.revised = case_when(
            cohort == "crossover" ~ cohort,
            subj.mixed.tracer ~ str_c(cohort, ".mixed_tracer"),
            !subj.mixed.tracer ~ str_c(cohort, ".single_tracer"),
        ),
        cohort.combined = str_c(cohort.revised, ".", tracer)
    )

# =============================
# ===== STATISTICAL TESTS =====
# =============================

# REVISED: tests performed:
# - head-to-head dataset:
#   - H2H vs. training set (all 4 cohorts)
#   - within training set
#       - FBP vs. PIB within single-tracer, mixed-tracer
#       - single-tracer vs. mixed-tracer within tracer
# - simulation dataset:
#   - FBP vs. PIB within single-tracer, mixed-tracer
#   - single-tracer vs. mixed-tracer within tracer

# ----- head-to-head dataset -----
head2head_df <- data_df %>%
    filter(cohort %in% c("crossover", "crossover_train"))

# test H2H vs. training set (all 4 cohorts)
h2h_train_vs_test <- head2head_df %>%
    filter(!(cohort == "crossover" & tracer == "AV45")) %>%  # remove one of the crossover cohorts (moved to single column for revised table)
    mutate(cohort = cohort.combined) %>%
    nest(.by = c(cohort.combined, tracer))
h2h_train_df <- h2h_train_vs_test %>%
    filter(cohort.combined == "crossover.PIB") %>%
    pull(data) %>% .[[1]]
h2h_train_vs_test <- h2h_train_vs_test %>%
    filter(cohort.combined != "crossover.PIB") %>%
    mutate(data = map(data, ~bind_rows(.x, h2h_train_df)))
h2h_train_vs_test_signif <- h2h_train_vs_test %>%
    rowwise() %>%
    mutate(
        test.age = data %>% t_test(age ~ cohort) %>% add_significance("p"),
        test.sex = table(data$cohort, data$sex) %>% fisher_test(),
        test.apoe = table(data$cohort, data$apoe) %>% fisher_test(),
        test.cdr = table(data$cohort, data$cdr) %>% fisher_test()
    ) %>%
    unnest(starts_with("test"), names_sep = "__") %>%
    select(cohort.combined, ends_with("p.signif")) %>%
    pivot_longer(
        cols = starts_with("test")
    ) %>%
    mutate(
        name = name %>% str_remove("test.") %>% str_remove("__p.signif")
    ) %>%
    pivot_wider(
        names_from = cohort.combined,
        values_from = value
    ) %>%
    rename_with(~ str_c(.x, ".trainvstest"), .cols = -name)

# test against tracers within cohort
head2head_tracer_test <- head2head_df %>%
    filter(cohort.revised != "crossover") %>%
    nest(.by = cohort.revised) %>%
    rowwise() %>%
    mutate(
        test.age = data %>% t_test(age ~ tracer) %>% add_significance("p"),
        test.sex = table(data$tracer, data$sex) %>% fisher_test(),
        test.apoe = table(data$tracer, data$apoe) %>% fisher_test(),
        test.cdr = table(data$tracer, data$cdr) %>% fisher_test()
    ) %>%
    unnest(starts_with("test"), names_sep = "__") %>%
    select(cohort.revised, ends_with("p.signif")) %>%
    pivot_longer(
        cols = starts_with("test")
    ) %>%
    mutate(
        name = name %>% str_remove("test.") %>% str_remove("__p.signif")
    ) %>%
    pivot_wider(
        names_from = cohort.revised,
        values_from = value
    ) %>%
    rename_with(~ str_c(.x, ".tracer_test"), .cols = -name) %>%
    mutate(across(
        -name,
        ~ str_replace_all(.x, "\\*", "\u271d")  # replace with cross
    ))

# ----- simulation dataset -----
simulation_df <- data_df %>%
    filter(cohort %in% c("simulation"))

simulation_test <- simulation_df %>%
    nest(.by = cohort.revised) %>%
    rowwise() %>%
    mutate(
        test.age = data %>% t_test(age ~ tracer) %>% add_significance("p"),
        test.sex = table(data$tracer, data$sex) %>% fisher_test(),
        test.apoe = table(data$tracer, data$apoe) %>% fisher_test(),
        test.cdr = table(data$tracer, data$cdr) %>% fisher_test()
    ) %>%
    unnest(starts_with("test"), names_sep = "__") %>%
    select(cohort.revised, ends_with("p.signif")) %>%
    pivot_longer(cols = starts_with("test")) %>%
    mutate(name = name %>% str_remove("test.") %>% str_remove("__p.signif")) %>%
    pivot_wider(
        names_from = cohort.revised,
        values_from = value
    ) %>%
    rename_with(~ str_c(.x, ".tracer_test"), .cols = -name) %>%
    mutate(across(
        -name,
        ~ str_replace_all(.x, "\\*", "\u271d")  # replace with cross
    ))

# ----- combine results -----
combine_signif <- function(signif1, signif2) {
    if_else(
        signif1 != "" & signif2 != "",
        str_c(signif1, signif2, sep = ","),
        str_c(signif1, signif2, sep = "")
    )
}
test_df <-
    reduce(list(h2h_train_vs_test_signif, head2head_tracer_test, simulation_test), ~ left_join(.x, .y, by = "name")) %>%
    mutate(
        across(
            -name,
            ~str_replace_all(.x, "ns", "")
        ),
        crossover_train.single_tracer.PIB.combined = combine_signif(crossover_train.single_tracer.PIB.trainvstest, crossover_train.single_tracer.tracer_test),
        crossover_train.mixed_tracer.PIB.combined = combine_signif(crossover_train.mixed_tracer.PIB.trainvstest, crossover_train.mixed_tracer.tracer_test), # this is for the "crossover training PiB" columns, since it will have results for both tests (across cohort, across tracer)
        name = if_else(name == "cdr", "n.clinical_group", name)
    )

# ============================
# ===== GET DEMOGRAPHICS =====
# ============================

# scan-level stats
# mean age at scan
# amyloid-positive scans
# closest assigned CDR
# mean time interval between scans
stats_scan <- data_df %>%
    group_by(tracer, cohort.revised) %>%
    summarise(
        across(
            c(age, age.diff),
            list("mean" = ~ mean(.x, na.rm = TRUE), "sd" = ~ sd(.x, na.rm = TRUE))
        ),
        n.scan = n(),
        n.amyloid_positive = sum(amyloid_positive),
        n.amyloid_negative = sum(!amyloid_positive),
        n.cn = sum(cdr == 0),
        n.mci = sum(cdr == 0.5),
        n.ad = sum(cdr >= 1),
        .groups = "drop"
    )

# subject-level stats
# number of subjects
# sex
# APOE
# mean number of scans per subject
stats_subj <- data_df %>%
    group_by(tracer, cohort.revised, subj) %>%
    summarise(
        n.scans.subj =  n(),
        sex = unique(sex),
        apoe = unique(apoe),
        .groups = "drop"
    ) %>%
    group_by(tracer, cohort.revised) %>%
    summarise(
        n.subj = n(),
        n.male = sum(sex == "M"),
        n.female = sum(sex == "F"),
        n.apoe1 = sum(apoe == 1, na.rm = TRUE),
        n.apoe0 = sum(apoe == 0, na.rm = TRUE),
        across(
            n.scans.subj,
            list("mean" = ~ mean(.x, na.rm = TRUE), "sd" = ~ sd(.x, na.rm = TRUE))
        ),
        .groups = "drop"
    )

# concatenate
demo_df <- 
    full_join(stats_scan, stats_subj, by = c("tracer", "cohort.revised")) %>%
    mutate(  # create string columns
        age = mean_sd_str(age_mean, age_sd),
        age.diff = mean_sd_str(age.diff_mean, age.diff_sd),
        sex = paste(n.male, n.female, sep = "/"),
        apoe = paste(n.apoe0, n.apoe1, sep = "/"),
        n.amyloid_positive = paste(n.amyloid_positive, n.amyloid_negative, sep = "/"),
        n.clinical_group = paste(n.cn, n.mci, n.ad, sep = "/"),
        n.scans.subj = mean_sd_str(n.scans.subj_mean, n.scans.subj_sd),
        .keep = "unused"
    ) %>%
    mutate(across(everything(), as.character))  # convert everything to characters

# reorder dataframe
f <- list(
    "n.subj" = "Number of subjects",
    "n.scan" = "Number of scans",
    "age" = "Mean age at scan (\u00B1 sd)",
    "sex" = "Number of males/females",
    "apoe" = "APOE noncarriers/carriers",
    "n.clinical_group" = "Number of CDR = 0/CDR = 0.5/CDR > 1",
    "n.scans.subj" = "Number of scans per subject (\u00B1 sd)",
    "age.diff" = "Mean years between scans (\u00B1 sd)"
)

append_significance <- function(c1, c2) {
    # append statistical test results
    if_else(
        c2 == "" | is.na(c2),
        c1,
        str_c(c1, "<sup>", c2, "</sup>")
    )
}

demo_df_reformat <- demo_df %>%
    select(cohort.revised, tracer, names(f)) %>%
    pivot_longer(
        cols = -c(tracer, cohort.revised)
    ) %>%
    pivot_wider(
        names_from = c(cohort.revised, tracer),
        values_from = value,
        names_sep = "."
    ) %>%
    left_join(test_df, by = "name") %>%
    mutate(
        crossover_train.single_tracer.AV45 = append_significance(crossover_train.single_tracer.AV45, crossover_train.single_tracer.AV45.trainvstest),
        crossover_train.single_tracer.PIB = append_significance(crossover_train.single_tracer.PIB, crossover_train.single_tracer.PIB.combined),
        crossover_train.mixed_tracer.AV45 = append_significance(crossover_train.mixed_tracer.AV45, crossover_train.mixed_tracer.AV45.trainvstest),
        crossover_train.mixed_tracer.PIB = append_significance(crossover_train.mixed_tracer.PIB, crossover_train.mixed_tracer.PIB.combined),
        simulation.single_tracer.PIB = append_significance(simulation.single_tracer.PIB, simulation.single_tracer.tracer_test),
        simulation.mixed_tracer.PIB = append_significance(simulation.mixed_tracer.PIB, simulation.mixed_tracer.tracer_test),
        across(ends_with(c("AV45", "PIB")), ~ map(.x, html))  # need to map each cell individually to HTML
    ) %>%
    mutate(name = f[name] %>% as.character()) %>%  # rename row names
    select(-ends_with("test"))

# REVISIONS: remove 1 column from crossover section
demo_df_reformat <- demo_df_reformat %>%
    select(-crossover.AV45) %>%
    rename(crossover = crossover.PIB) %>%
    select(name, crossover,
           crossover_train.single_tracer.AV45, crossover_train.single_tracer.PIB,
           crossover_train.mixed_tracer.AV45, crossover_train.mixed_tracer.PIB,
           simulation.single_tracer.AV45, simulation.single_tracer.PIB,
           simulation.mixed_tracer.AV45, simulation.mixed_tracer.PIB,
    )

# =====================================
# ===== CREATE DEMOGRAPHICS TABLE =====
# =====================================

demographics_gt <- demo_df_reformat %>%
    gt(rowname_col = "name") %>%
    cols_label(
        crossover = "",
        crossover_train.single_tracer.AV45 = "FBP", crossover_train.mixed_tracer.AV45 = "FBP",
        simulation.single_tracer.AV45 = "FBP", simulation.mixed_tracer.AV45 = "FBP",
        crossover_train.single_tracer.PIB = "PiB", crossover_train.mixed_tracer.PIB = "PiB",
        simulation.single_tracer.PIB = "PiB", simulation.mixed_tracer.PIB = "PiB"
    ) %>%
    tab_spanner(
        label = "head-to-head",
        columns = c("crossover"),
        level = 2
    ) %>%
    tab_spanner(
        label = "single-tracer",
        columns = c("crossover_train.single_tracer.AV45", "crossover_train.single_tracer.PIB"),
        id = "single_tracer_h2h"
    ) %>%
    tab_spanner(
        label = "mixed-tracer",
        columns = c("crossover_train.mixed_tracer.AV45", "crossover_train.mixed_tracer.PIB"),
        id = "mixed_tracer_h2h"
    ) %>%
    tab_spanner(
        label = "training",
        spanners = c("single_tracer_h2h", "mixed_tracer_h2h")
    ) %>%
    tab_spanner(
        label = "Tracer head-to-head comparison",
        spanners = c("head-to-head", "training")
    ) %>%
    tab_spanner(
        label = "single-tracer",
        columns = c("simulation.single_tracer.AV45", "simulation.single_tracer.PIB"),
        id = "single_tracer_sim"
    ) %>%
    tab_spanner(
        label = "mixed-tracer",
        columns = c("simulation.mixed_tracer.AV45", "simulation.mixed_tracer.PIB"),
        id = "mixed_tracer_sim"
    ) %>%
    tab_spanner(
        label = "Simulation",
        spanners = c("single_tracer_sim", "mixed_tracer_sim"),
        level = 3
    ) %>%
    tab_spanner(
        label = "",
        spanners = c("Simulation")
    ) %>%
    tab_footnote(
        footnote = "sd = standard deviation",
        location = cells_stub(
            rows = "Mean age at scan (\u00B1 sd)"
        )
    ) %>%
    tab_footnote(
        footnote = "APOE = apolipoprotein-E4",
        location = cells_stub(
            rows = "APOE noncarriers/carriers"
        )
    ) %>%
    tab_footnote(
        footnote = "CDR = Clinical Dementia Rating",
        location = cells_stub(
            rows = "Number of CDR = 0/CDR = 0.5/CDR > 1"
        )
    ) %>%
    tab_footnote("* = significance of head-to-head vs. training comparison") %>%
    tab_footnote("\u271d = significance of tracer vs. tracer comparison")

# ==============================
# ===== SAVE TABLE TO FILE =====
# ==============================

gtsave(demographics_gt, "demographics.html", opt$odir)
gtsave(demographics_gt, "demographics.rtf", opt$odir)

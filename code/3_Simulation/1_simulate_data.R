# ==============================================================================

# Name: 1_simulate_data.R
# Author: Braden Yang
# Created: 12/20/2023
# Description: Simulate longitudinal amyloid PET data with LMER models trained on
#   OASIS data, and introduce engineered treatment effect

# Notes
# - simulated data is meant to mimic the data acquired from a successful clinical
#   trial, where an observable difference in amyloid rate-of-change exists
#   between "treatment" and "placebo" groups
# - both the "treatment" and "placebo" group data will be generated from an LME
#   model learned on cognitively normal, amyloid-positive subjects from OASIS. The
#   only difference is that the "treatment" group will receive an additional
#   treatment effect term, which will be determined from Charlie's annualized
#   rate-of-change data for each region
# - step-by-step procedure
#   1. select cohort from OASIS to train LME
#   2. extract all features (or select features of interest)
#   3. pivot table, nest by feature and tracer
#   4. fit LME for each feature and tracer
#   5. generate data for 2 groups, use longitudinal treatment effect for treatment group
#   6. pivot wider, save as CSV

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
    make_option(c("-o", "--odir"), action="store", default="data/simulation",
        type="character", help="Path to output directory"),
    make_option(c("-n", "--n_simulation"), action="store", type="integer", default=1000,
        help="number of simulations to perform"),
    make_option(c("-p", "--p_av45_1"), action="store", type="double", default=0.5,
        help="proportion of AV45 scans in bootstrap, given in decimal fraction"),
    make_option(c("-q", "--p_av45_2"), action="store", type="double", default=NULL,
        help="(optional) proportion of AV45 for 2nd clinical group. Specify if you want different tracer proportions across clinical groups. If not given, the proportion will be made the same across groups, using p_av45_1 for both"),
    make_option(c("-s", "--n_subj"), action="store", type="integer", default=30,
        help="number of subjects per simulated group"),
    make_option(c("-t", "--treatment.delta"), action="store", type="double", default=0.01,
        help="annualized rate-of-change of the treatment group"),
    make_option(c("-r", "--random_seed"), action="store", type="integer", default=NULL,
        help="random seed for bootstrap sampling"),
    make_option(c("-b", "--assign_tracer_by_subj"), action="store_true", default=FALSE,
        help="if specified, do subject-wise random tracer assignment rather than scan-wise; p_av45 will refer to the proportion of subjects assigned AV45 in this case"),
    make_option(c("-c", "--pvc"), action="store_true", default=FALSE,
        help="if specified, run simulations on PVC data"),
    make_option(c("-a", "--save_amypos"), action="store_true", default=FALSE,
        help="if specified, just save amypos.csv (for training ComBat) without running simulations"),
    make_option(c("-x", "--alt_tracer_assign"), action="store_true", default=FALSE,
        help="use a more realistic method for assigning tracers, e.g. have the subject switch tracers only once"),
    make_option(c("-u", "--var_multiplier_int"), action="store", type="double", default=1,
        help="multiplier for variance of random intercept"),
    make_option(c("-v", "--var_multiplier_res"), action="store", type="double", default=1,
        help="multiplier for variance of residual")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (INTERACTIVE) {
    opt$wdir <- "/home/b.y.yang/sotiraslab/BradenPETHarmonization"
    opt$n_simulation <- 10
    opt$treatment.delta <- 0.03
    opt$n_subj <- 50
    opt$p_av45_1 <- 0.2
    opt$p_av45_2 <- 0.6
    opt$random_seed <- 42
    opt$assign_tracer_by_subj <- FALSE
    opt$pvc <- FALSE
}

# set working directory
if (!is.null(opt$wdir)) {
    setwd(opt$wdir)
}

# load byyfunctions
load_all("submodules/byyfunctions")

# make output directory
byyfunctions::make_dir(opt$odir)

# ======================================
# ===== SOURCE FUNCTIONS/VARIABLES =====
# ======================================

source("code/functions.R")

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

get_time_from_baseline <- function(.data, age_col, group_col = "subj") {

    df <- .data %>%
        group_by(pick(all_of(group_col))) %>%
        mutate(
            age.baseline = min(.data[[age_col]]),
            time.from.baseline = .data[[age_col]] - age.baseline
        ) %>%
        ungroup()

    return(df)

}

simulate_trial <- function(n_subj, p_av45_1, p_av45_2 = NULL, random_seed = NULL, assign_tracer_by_subj = FALSE, realistic = FALSE, realistic_start_tracer = NULL) {

    # helper functions
    assign_tracer <- function(.data, p, group_col = NULL, time_col = null, assign_by_subj = FALSE, realistic = FALSE, realistic_start_tracer = NULL) {

        # randomly assign tracer to each row, according to
        # the specified proportion `p`
        # 
        # if `assign_by_subj=TRUE`, then the dataframe will
        # first be nested by `group_col`, then the assignment
        # will be performed on subjects. Thus, `p` will equate
        # to the proportion of *subjects* who are assigned AV45
        # 
        # if `realistic=TRUE`, then tracers will be assigned
        # in a realistic way, e.g. subjects will switch tracers
        # only once

        sample_prop <- function(x, n, p = NULL, n_sample = NULL) {

            # generate binomial permutation with (almost) exact proportion of trues to falses

            # if n = 0, return null
            if (n == 0) {
                return(NULL)
            }

            if (is.null(n_sample)) {
                n_sample <- round(n*p)
            }

            idx <- sample(1:n, size = n_sample, replace = FALSE)

            v <- if_else(
                1:n %in% idx,
                x[[1]],
                x[[2]]
            )

            return(v)

        }

        order_tracers <- function(x, start_tracer) {
            start_tracer <- start_tracer[[1]]
            if (start_tracer == "AV45") {
                l <- c("AV45", "PIB")
            } else {
                l <- c("PIB", "AV45")
            }
            x_sort <- x %>% ordered(levels = l) %>% sort %>% as.character
            return(x_sort)
        }

        if (assign_by_subj) {
            .data <- .data %>%
                nest(.by = all_of(group_col)) %>%
                mutate(tracer = sample_prop(c("AV45", "PIB"), n(), p)) %>%
                unnest(data)
        } else {
            .data <- .data %>%
                mutate(tracer = sample_prop(c("AV45", "PIB"), n(), p))
            if (realistic) {
                .data <- .data %>%
                    group_by(pick(all_of(group_col))) %>%
                    arrange(.data[[time_col]]) %>%
                    ungroup() %>%
                    nest(.by = all_of(group_col)) %>%
                    mutate(
                        .start_tracer = if (is.null(realistic_start_tracer)) {sample(c("AV45", "PIB"), size = n(), replace = TRUE)} else {realistic_start_tracer},
                    ) %>%
                    unnest(data) %>%
                    group_by(pick(all_of(group_col))) %>%
                    mutate(
                        tracer = order_tracers(tracer, .start_tracer)
                    ) %>%
                    ungroup() %>%
                    select(-.start_tracer)
            }
        }
        
        return(.data)

    }

    sample_check <- function(x, size, replace = FALSE) {

        s <- rep(0, size)
        while(length(unique(s)) <= 1) {
            s <- sample(x, size, replace = replace)
        }
        return(s)

    }

    generate_timepoints <- function(n_timepoints_v, baseline_age_v, interval_v, random_seed = NULL) {

        # generate a random list of ages, based on the empirical distribution of baseline ages,
        # empirical list of possible time between scan, and empirical list of # of timepoints

        # set random seed
        if (!is.null(random_seed)) {
            set.seed(random_seed)
        }

        n_timepoints <- sample(n_timepoints_v, 1)
        baseline_age <- sample(baseline_age_v, 1)
        time_interval <- c(0, sample(interval_v, n_timepoints - 1, replace = TRUE))
        timepoints <- baseline_age + cumsum(time_interval)

        return(timepoints)

    }

    # set random seed
    if (!is.null(random_seed)) {
        set.seed(random_seed)
    }

    # check if p_av45_2 is NULL
    if (is.null(p_av45_2)) {
        p_av45_2 <- p_av45_1
    }

    # simlulate subject timepoints
    sim_data <- tibble(subj = str_c("sim_subj_", 1:(2*n_subj))) %>%
        mutate(
            trial.group = c(rep(0, n_subj), rep(1, n_subj)),
            trial.group.fct = factor(trial.group, levels = c(0, 1), labels = c("placebo", "treatment")),
            random_seed.timepoint = sample(c(1:1000000), n(), replace = FALSE),
            age = map(
                random_seed.timepoint,
                ~ generate_timepoints(scans_per_subj, age.baseline_v, time_between_scans, .x)
            ),
            sex = sample_check(sex_v, n(), replace = TRUE),
            apoe = sample_check(apoe_v, n(), replace = TRUE)
        ) %>%
        unnest(age) %>%
        nest(.by = trial.group) %>%
        left_join(tibble(trial.group = c(0, 1), p_av45 = c(p_av45_1, p_av45_2)), by = "trial.group") %>%
        rowwise() %>%
        mutate(
            data = data %>%
                assign_tracer(p_av45, "subj", "age", assign_tracer_by_subj, realistic, realistic_start_tracer) %>%
                list()
        ) %>%
        unnest(data) %>%
        get_time_from_baseline(age_col = "age", group_col = "subj") %>%
        mutate(time.from.baseline_trial.group = time.from.baseline * trial.group) %>%
        nest(.by = "tracer", .key = "sample")

    # generate SUVR data from LMER models
    lmer_sim_data <- lmer_model_df %>%
        left_join(sim_data, by = "tracer") %>%
        rowwise() %>%
        mutate(
            sample = bind_cols(
                sample,
                simulate(
                    lmer_model,
                    nsim = 1,
                    newdata = sample,
                    re.form = NA,
                    allow.new.levels = TRUE
                ) %>% as_tibble()
            ) %>% list()
        ) %>%
        select(tracer, feature_roi, sample) %>%
        unnest(sample) %>%
        rename(!!resp_var := sim_1)

    # add treatment effect
    sim_trial_data <- lmer_sim_data %>%
        bind_cols(delta_df) %>%
        mutate(
            delta = case_match(trial.group, 0 ~ placebo.delta, 1 ~ treatment.delta, .default = 0),
            trial.effect = delta * time.from.baseline,
            feature_value = feature_value + trial.effect
        )

    return(sim_trial_data)

}

# ============================
# ===== DEFINE VARIABLES =====
# ============================

# linear mixed effects formulas
resp_var <- "feature_value"
fixed_eff <- c("1", "age.baseline", "time.from.baseline", "sex", "apoe")
random_eff <- "(1 | subj)"
lhs <- paste(c(fixed_eff, random_eff), collapse = " + ")
lmer_formula <- paste(c(resp_var, "~", lhs), collapse = " ")

# load feature names
feature_names <- read_csv("data/csv/feature_names_trunc.csv", col_names = FALSE, show_col_types = FALSE) %>% pull(1)
feature_names <- feature_names[!str_detect(feature_names, "cerebellum_cortex")]  # remove cerebellum cortex

# set random seed
if (!is.null(opt$random_seed)) {
    set.seed(opt$random_seed)
}

# get random seeds for simulation
random_seed_v <- sample(c(1:1000000), opt$n_simulation, replace = FALSE)

# =================================
# ===== DEFINE RATE-OF-CHANGE =====
# =================================

# keep the rate-of-change constant across regions

delta_df <- tibble(
    tibble(
        placebo.delta = 0,
        treatment.delta = -1 * opt$treatment.delta
    )
)

# =====================
# ===== LOAD DATA =====
# =====================

if (opt$pvc) {
    # suvr_path <- "data/baseline_pvc.RDS"
} else {
    suvr_path <- "data/pet.RDS"
}

# load SUVR data
pet_df <- read_rds(suvr_path) %>%
    mutate(
        harmonization_method = "suvr",
        age = round_age(age)
    ) %>%
    filter(
        study == "OASIS",
        pvc == opt$pvc
    )

# =====================================
# ===== GET SIMULATION PARAMETERS =====
# =====================================

# get distributions of # of scans per subject, time in between scans and age at baseline
scans_filter <- pet_df %>%
    filter(harmonization_method == "suvr") %>%
    byyfunctions::get_diff_adjacent("age", "subj") %>%
    filter(is.na(age.diff) | age.diff > 0.5)  # remove scans that occur in close proximity
scans_per_subj <- scans_filter %>%
    count(subj) %>%
    filter(n > 1) %>%
    pull(n)
time_between_scans <- scans_filter %>%
    drop_na(age.diff) %>%
    pull(age.diff)
subj_stats <- pet_df %>%
    group_by(subj) %>%
    summarise(
        age.baseline = min(age, na.rm = TRUE),
        sex = unique(sex),
        apoe = unique(apoe)
    )
age.baseline_v <- subj_stats %>% pull(age.baseline)
sex_v <- subj_stats %>% pull(sex)
apoe_v <- subj_stats %>% drop_na(apoe) %>% pull(apoe)

# ================================
# ===== FILTER AND TIDY DATA =====
# ================================

select_pet <- function(.data) {
    df <- .data %>%
        select(
            subj, age, tracer, harmonization_method, sex, apoe, cdr, amyloid_positive,
            all_of(feature_names)
        )
    return(df)
}

# amyloid-positive scans
amypos_df <- pet_df %>%
    group_by(subj, tracer) %>%
    filter(any(amyloid_positive)) %>%
    ungroup() %>%
    select_pet() %>%
    drop_na(age, sex, apoe)

# save CSV of amyloid-positive data to train ComBat
if (opt$save_amypos) {

    write_csv(
        amypos_df %>% get_time_from_baseline(age_col = "age", group_col = "subj"),
        file.path(opt$odir, "amypos.csv")
    )

    stop("++ amypos.csv saved; exiting ++")
    
}

# pivot to long format
long_df <- amypos_df %>%
    pivot_longer(
        cols = all_of(feature_names),
        names_to = "feature_roi",
        values_to = "feature_value"
    )

# get other parameters
long_df <- long_df %>%
    get_time_from_baseline(age_col = "age", group_col = "subj")  # get time from baseline, baseline age

# ====================
# ===== FIT LMER =====
# ====================

lmer_model_df <- long_df %>%
    nest(.by = c(tracer, feature_roi)) %>%
    mutate(
        lmer_model = map(
            data,
            ~ lmer(
                formula = formula(lmer_formula),
                data = .x
            )
        )
    )

# =========================
# ===== SIMULATE DATA =====
# =========================

sim_trial_df <-
    tibble(
        sim_idx = 1:length(random_seed_v),
        random_seed.sim = random_seed_v
    ) %>%
    mutate(sim_trial_stats = map(
        random_seed.sim,
        ~ simulate_trial(
            n_subj = opt$n_subj,
            p_av45_1 = opt$p_av45_1,
            p_av45_2 = opt$p_av45_2,
            random_seed = .x,
            assign_tracer_by_subj = opt$assign_tracer_by_subj,
            realistic = opt$alt_tracer_assign
        )
    )) %>%
    unnest(sim_trial_stats) %>%
    pivot_wider(
        names_from = feature_roi,
        values_from = feature_value
    )

# save results
write_csv(
    sim_trial_df,
    file.path(opt$odir, "sim.csv")
)

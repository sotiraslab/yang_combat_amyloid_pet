# ==============================================================================

# Name: get_data.R
# Author: Braden Yang
# Created: 06/05/2023
# Description: get cross-sectional amyloid PET features, tidy, and match clinical
#   and other non-imaging data

# ==============================================================================

# clear environment
rm(list = ls())

# *** toggle variables ***
INTERACTIVE <- FALSE
# *** toggle variables ***

# ===========================
# ===== IMPORT PACKAGES =====
# ===========================

library(rjson)
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
devtools::load_all("submodules/byyfunctions")

# ============================
# ===== DEFINE FUNCTIONS =====
# ============================

oasis_metadata_from_id <- function(.data, id_col, sep = "_", num_sep = 3) {

    # Return dataframe with new columns "Subject" and "daysFromEntry",
    # infered from the column "PUP_PUPTIMECOURSEDATA.ID"
    # 
    # Parameters
    # ----------
    #   .data:   OASIS-3 dataframe of PUP data (loaded from `oasis3_pup_full.csv`)
    #   id_col:  column containing unique subject/session ID (use tidyverse data
    #            masking when inputing column, i.e. no quotes)
    #   sep:     separating character (default = "_")
    #   num_sep: how many substrings result after splitting by sep?
    # 
    # Returns
    # -------
    #   dataframe with new columns

    # split string into components
    id_split <- .data %>%
        pull({{id_col}}) %>%
        str_split_fixed(pattern = sep, n = num_sep)

    # get subject ID as new column "Subject"
    .data$Subject <- id_split[, 1]

    # get days from entry as new column "daysFromEntry"
    # remove first character of every string: https://stackoverflow.com/questions/35113553/r-remove-first-character-from-string
    .data$daysFromEntry <- as.integer(sub(".", "", id_split[, num_sep]))

    return(.data)

}

oasis_match_features <- function(pet_df, cdr_df, fs_df, demog_df) {

    # ===== 1. Get APOE, age and sex =====

    pet_df <- demog_df %>%
        rename(Subject = OASISID, apoe = APOE, ageAtEntry = AgeatEntry, sex = GENDER) %>%
        mutate(
            sex = factor(sex, levels = c(1, 2), labels = c("M", "F")),
            apoe = as.factor(
                case_when(
                    apoe %in% c("34", "44") ~ 1,
                    apoe %in% c("22", "23", "24", "33") ~ 0,
                    .default = NA
                )
            ),
        ) %>%
        group_by(Subject) %>%
        slice_head(n = 1) %>%
        select(Subject, apoe, ageAtEntry, sex) %>%
        right_join(pet_df, by = "Subject", multiple = "all")

    # ===== 2. Merge FreeSurfer volumetric data =====

    pet_df <- fs_df %>%  # then merge
        select(FSId, freesurfer_version, mri_field_strength, all_of(roi_features_df %>% filter(study == "oasis", feature_type == "vol") %>% pull(col))) %>%
        right_join(pet_df, by = "FSId", multiple = "all")

    # ===== 3. Compute age at visit =====

    pet_df <- pet_df %>%
        mutate(ageAtVisit = ageAtEntry + (daysFromEntry) / (365.25))

    # ===== 4. Match CDR and MMSE =====

    pet_df <- pet_df %>%
        byyfunctions::match_data(
            match_df = cdr_df,
            ref_col = "daysFromEntry",
            match_col = "daysFromEntry",
            select_col = c("cdr", "mmse"),
            group_col = "Subject"
        ) %>%
        rename(
            cdr = daysFromEntry_match_cdr,
            mmse = daysFromEntry_match_mmse
        )

    # ===== 5. Rename columns =====

    pet_df <- pet_df %>%
        byyfunctions::rename_columns(
            col_old = roi_features_df %>% filter(study == "oasis") %>% pull(col),
            col_new = roi_features_df %>% filter(study == "oasis") %>% pull(fs_label)
        ) %>%
        rename(
            subj = Subject,
            age = ageAtVisit,
            summary_oasis.suvr = PET_fSUVR_TOT_CORTMEAN
        ) %>%
        select(-c(ageAtEntry, FSId), -starts_with("daysFromEntry"))

    return(pet_df)

}

mark_amyloid_positive <- function(suvr, tracer, study, pvc) {

    # for ADNI, use the summary SUVR w/whole cerebellum normalization
    # for OASIS, use mean cortical SUVR

    amyloid_positive <- case_when(
        study == "ADNI" & tracer == "AV45" ~ suvr >= 1.11,
        study == "ADNI" & tracer == "FBB" ~ suvr >= 1.08,
        study == "OASIS" & tracer == "AV45" & !pvc ~ suvr >= 1.24,
        study == "OASIS" & tracer == "PIB" & !pvc ~ suvr >= 1.31,
        study == "OASIS" & tracer == "AV45" & pvc ~ suvr >= 1.19,
        study == "OASIS" & tracer == "PIB" & pvc ~ suvr >= 1.42,
        .default = NA
    )

    return(amyloid_positive)

}

assign_clinical_group <- function(cdr_vec) {

    group_vec <- dplyr::case_when(
        {{cdr_vec}} == 0 ~ "CN",
        {{cdr_vec}} == 0.5 ~ "MCI",
        {{cdr_vec}} > 0.5 ~ "AD",
        .default = NA 
    ) %>%
    factor(levels = c("CN", "MCI", "AD"))

    return(group_vec)

}

extend_clinical_group <- function(clinical_group, amyloid_positive) {

    group_vec <- dplyr::if_else(
        {{clinical_group}} == "CN" & {{amyloid_positive}},
        "preclinical",
        {{clinical_group}}
    ) %>%
    factor(levels = c("CN", "preclinical", "MCI", "AD"))

    return(group_vec)

}

normalize_by_reference <- function(.data, cols_to_norm, ref_cols) {

    # Divide a set of columns (cols_to_norm) by the average of the list of
    # columns defined in ref_cols, return the normalized table

    # summing multiple columns: https://stackoverflow.com/questions/47759347/create-a-new-column-which-is-the-sum-of-specific-columns-selected-by-their-name

    .data <- .data %>%
        mutate(ref = reduce(select(., all_of(ref_cols)), `+`) / length(ref_cols)) %>%
        mutate_at(cols_to_norm, ~ . / ref) %>%
        select(-ref)
    
    return(.data)

}

assigned_grouped_idx <- function(.data, .by) {

    # create row index column for a grouped tibble
    # new column is named `idx`

    .data <- .data %>%
        nest(.by = {{.by}}) %>%
        mutate(idx = row_number()) %>%
        unnest(data)

    return(.data)

}

get_baseline <- function(.data, group_col, time_col) {

    .data <- .data %>%
        group_by({{group_col}}) %>%
        slice_min({{time_col}}) %>%
        ungroup()

    return(.data)

}

# =======================
# ===== LOAD PARAMS =====
# =======================

params <- fromJSON(file = "params.json")
oasis_paths <- params$OASIS_paths

# =======================
# ===== LOAD TABLES =====
# =======================

print("++ Loading tables ++")

# load tables
oasis_df <- purrr::map(oasis_paths, ~read_csv(.x, show_col_types = FALSE))

# table of ROI features
roi_features_df <- 
    read_rds("data/csv/roi.RDS") %>%
    mutate(col = case_when(
        study == "adni" ~ str_remove(col, ".av45"),
        study == "oasis" & feature_type == "suvr" ~ str_remove(col, ".pup"),
        study == "oasis" & feature_type == "vol" ~ str_replace_all(str_remove(col, ".fs"), "\\.", "-")
    )) %>%
    mutate(fs_label = paste(fs_label, feature_type, sep = "."))

# get names of features
feature_names <- read_csv("data/csv/feature_names.csv", col_names = FALSE) %>% pull(1)

# ==================================
# ===== EXTRACT OASIS FEATURES =====
# ==================================

print("++ Preparing OASIS data ++")

# get new columns (Subject, daysFromEntry) in clinical_df and pup_full_df
oasis_pet_df <- oasis_df[["pet"]] %>%
    oasis_metadata_from_id(id_col = `PUP_PUPTIMECOURSEDATA ID`, num_sep = 4)
oasis_cdr_df <- oasis_df[["cdr"]] %>%
    rename(Subject = OASISID, daysFromEntry = days_to_visit, cdr = CDRTOT, mmse = MMSE)

# REVISIONS: get list of WM regions
wm_roi <- oasis_df[["pet"]] %>%
    select(starts_with(c("PET_fSUVR_L_WM_", "PET_fSUVR_R_WM_")), -contains(c("UNSEGMENTED", "CRPCLM"))) %>%
    colnames

# REVISIONS: get MRI df for magnetic field strength info
oasis_mri_df <- oasis_df[["mri_meta"]] %>%
    filter(`scan category` == "T1w") %>%
    group_by(label) %>% slice(1) %>% ungroup() %>%  # remove duplicate rows (assumed to be different runs)
    select(label, MagneticFieldStrength) %>%
    rename(MR_session = label, mri_field_strength = MagneticFieldStrength) %>%
    mutate(mri_field_strength = round(mri_field_strength, 1))  # some 1.49T exist for some reason...

# REVISIONS: get PET scanner model
oasis_pet_meta_df <- oasis_df[["pet_meta"]] %>%
    select(session_id, ManufacturersModelName) %>%
    rename(pet_session_id = session_id, pet_scanner_model = ManufacturersModelName) %>%
    mutate(pet_scanner_model = if_else(pet_scanner_model == "Biograph 40 PET CT", "Biograph 40 PET/CT", pet_scanner_model))

# get new columns in FS df and rename columns
oasis_fs_df <- oasis_df[["freesurfer"]] %>%
    rename(FSId = `FS_FSDATA ID`) %>%
    mutate(
        daysFromEntry = stringr::str_split_fixed(FSId, pattern = "_", n=3)[,3] %>% sub(".", "", .) %>% as.integer(),
        freesurfer_version = stringr::str_extract(version, "v\\d\\.\\d\\.\\d")
    ) %>%
    left_join(oasis_mri_df, by = "MR_session")

# select relevant cols in oasis_pet_df
oasis_suvr_features <- roi_features_df %>%
    filter(study == "oasis", feature_type == "suvr") %>%
    pull(col)
if (opt$white_matter) {
    oasis_suvr_features <- c(oasis_suvr_features, wm_roi)
}
oasis_pvc_suvr_features <- oasis_suvr_features %>%
    str_replace("PET_fSUVR", "PET_fSUVR_rsf") %>%
    str_remove(".pup")
oasis_pet_df <- oasis_pet_df %>%
    mutate(pet_session_id = str_replace(`PUP_PUPTIMECOURSEDATA ID`, "_PUPTIMECOURSE", "")) %>%  # REVISIONS: merge PET scanner model
    left_join(oasis_pet_meta_df, by = "pet_session_id") %>%
    select(
        Subject, daysFromEntry, FSId, tracer, pet_scanner_model,
        all_of(c(oasis_suvr_features, oasis_pvc_suvr_features)),
        PET_fSUVR_TOT_CORTMEAN, PET_fSUVR_rsf_TOT_CORTMEAN
    )

# pivot longer for PVC/non-PVC
oasis_pet_df <- oasis_pet_df %>%
    pivot_longer(
        cols = starts_with("PET_fSUVR"),
        names_pattern = "PET_fSUVR_(rsf_|)(.*)",
        names_to = c("pvc", "feature_roi"),
        values_to = "feature_value"
    ) %>%
    mutate(feature_roi = paste0("PET_fSUVR_", feature_roi)) %>%
    pivot_wider(
        names_from = "feature_roi",
        values_from = "feature_value"
    ) %>%
    mutate(pvc = ifelse(pvc == "rsf_", TRUE, FALSE))

# match features
oasis_pet_df <- oasis_pet_df %>%
    oasis_match_features(
        cdr_df = oasis_cdr_df,
        fs_df = oasis_fs_df,
        demog_df = oasis_df[["demog"]]
    )

# ====================================
# ===== ADDITIONAL PREPROCESSING =====
# ====================================

oasis_pet_df <- oasis_pet_df %>%
    mutate(study = "OASIS") %>%
    filter(if_all(c(age, sex, cdr), ~!is.na(.x))) %>%
    rename(summary.suvr = summary_oasis.suvr) %>%
    mutate(
        amyloid_positive = mark_amyloid_positive(summary.suvr, tracer, study, pvc),  # mark amyloid-positivity by summary cortical region
        clinical_group = assign_clinical_group(cdr),  # mark clinical groups by CDR
        clinical_group_extended = extend_clinical_group(clinical_group, amyloid_positive)  # mark preclinical subjects
    ) %>%
    select(-ends_with(c(".vol", "_oasis.suvr"))) %>%  # remove and rename columns
    rename_with(~str_remove(.x, ".suvr"), ends_with(".suvr")) %>%  # remove suffix from colnames
    filter(if_all(all_of(c(feature_names, "summary")), ~!is.na(.x))) %>%  # remove any rows with NA in specified ROIs
    assigned_grouped_idx(c(subj, age)) %>%  # assign unique index to each row (group by subj and age)
    normalize_by_reference(  # normalize by reference region 
        feature_names[! feature_names %in% "summary"],  # don't normalize summary; this is already normalized for us
        c("left.cerebellum.cortex", "right.cerebellum.cortex")  # average of left and right cerebellar cortex as ref
    )

# get cross-sectional baseline sampling
baseline_df <- oasis_pet_df %>% filter(!pvc) %>% get_baseline(subj, age)
baseline_pvc_df <- oasis_pet_df %>% filter(pvc) %>% get_baseline(subj, age)

# =======================
# ===== SAVE TABLES =====
# =======================

print("++ Saving RDS file ++")

if (opt$white_matter) {
    out_suffix <- "_wm"
} else {
    out_suffix <- ""
}

# create output directory
byyfunctions::make_dir(opt$odir, verbose = FALSE)

# write pet_df to RDS
write_rds(oasis_pet_df, file.path(opt$odir, str_glue("pet{out_suffix}.RDS")))
write_rds(baseline_df, file.path(opt$odir, str_glue("baseline{out_suffix}.RDS")))
write_rds(baseline_pvc_df, file.path(opt$odir, str_glue("baseline_pvc{out_suffix}.RDS")))

# REVISIONS: get WM regions and save as CSV
if (opt$white_matter && !file.exists(file.path(opt$odir, "csv/white_matter_regions.csv"))) {
    byyfunctions::make_dir(file.path(opt$odir, "csv"))
    wm_roi %>%
        as_tibble %>%
        write_csv(file.path(opt$odir, "csv/white_matter_regions.csv"), col_names=FALSE)
}

#' @title Generate roi_features.RData for byyfunctions package
#' 
#' @description Internal function used to generate the dataset `roi_features.RData`,
#' which contains a tibble of FreeSurfer cortical/subcortical ROI from the DK atlas,
#' along with their corresponding rows in the ADNI/OASIS merged tables and ggseg-
#' compatible labels
#' 
#' @param roi_features_csv path to ROI features CSV file
#' @param outfile path of output RData file containing tibble
#' 
#' @noRd
create_roi_features <- function(

    roi_features_csv = "/home/b.y.yang/BradenADLongitudinalPrediction/data/adni_oasis_roi_features.csv",
    summary_cort_roi_csv = "/home/b.y.yang/BradenADLongitudinalPrediction/data/summary_cortical_roi.csv",
    outfile = "/home/b.y.yang/BradenADLongitudinalPrediction/code/byyfunctions/data/roi_features.RData"
    ) {

    # load table of ROI feature names
    roi_features_raw <- readr::read_csv(roi_features_csv, col_types="ccccc") %>%
        dplyr::mutate(
            adni_suvr = byyfunctions::paste0_na(adni_suvr, ".av45"),
            oasis_suvr = byyfunctions::paste0_na(oasis_suvr, ".pup"),
            adni_vol = byyfunctions::paste0_na(adni_vol, ".av45"),
            oasis_vol = byyfunctions::paste0_na(oasis_vol, ".fs")
        ) # append suffix to each column label

    # tidy data
    roi_features <- roi_features_raw %>%
        tidyr::pivot_longer(
            cols = -"fs_label",
            names_pattern = "(.*)_(.*)",
            names_to = c("study", "feature_type"),
            values_to = "col"
        ) %>%
        tidyr::drop_na() %>%  # remove rows with NA
        byyfunctions::add_ggseg_label("fs_label")  # map ggseg labels to each ROI feature name

    # categorize each ROI (ADNI summary cortical, other cortical, subcortical)
    summary_roi_df <- readr::read_csv(summary_cort_roi_csv) %>%
        dplyr::mutate(fs_label = roi)
    roi_type_labels <- c("Summary cortical ROI", "Other cortical ROI", "Subcortical ROI")
    roi_features <- roi_features %>%
        dplyr::mutate(
            roi_type = ifelse(
                fs_label %in% summary_roi_df$fs_label,
                1,
                ifelse(startsWith(fs_label, "ctx."), 2, 3))
        ) %>%
        dplyr::mutate(roi_type = factor(roi_type, levels = c(1, 2, 3), labels = roi_type_labels))
    
    # if summary cortical, add lobe
    roi_features <- roi_features %>%
        dplyr::left_join(summary_roi_df %>% dplyr::select(fs_label, lobe), by = "fs_label")

    # save in "code/byyfunctions/data"
    save(roi_features, file = outfile)

}

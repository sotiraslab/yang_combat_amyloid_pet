#' @title Select and rename FreeSurfer ROI-based features from ADNI or OASIS merged table
#' 
#' @description This function first selects columns corresponding to FreeSurfer 
#' ROI-based features from the ADNI or OASIS datasets, along with other user-
#' specified columns, then renames the FS columns to FreeSurfer standard ROI
#' labels as notated in the color lookup table
#
#' @param .data ADNI or OASIS merged dataframe
#' @param col_select character vector containing additional columns to select
#' @param study_str indicates study, either "adni" or "oasis" (must be lowercase);
#' defaults to adni
#' @param feature_type_str (optional) indicates feature type, either "suvr" or "vol";
#' if not specified, then both feature types will be selected
#
#' @return dataframe with selected and renamed columns
#' 
#' @export
select_and_rename_fs <- function(.data, col_select, study_str = "adni", feature_type_str = NULL) {

    # get study and feature-type specific column names
    filter_df <- byyfunctions::roi_features %>%
        dplyr::filter(
            study == study_str,
            if (is.null(feature_type_str)) {TRUE} else {feature_type == feature_type_str}
        ) %>%
        dplyr::mutate(fs_label = paste(fs_label, feature_type, sep = ".")) # append feature_type to fs_label
    names_from <- filter_df$col
    names_to <- filter_df$fs_label
    
    .data <- .data %>%
        dplyr::select(dplyr::all_of(c(col_select, names_from))) %>% # first select columns
        byyfunctions::rename_columns(names_from, names_to) # then rename

    return(.data)

}
#' @title Add new column with ggseg ROI labels
#' 
#' @description This function transforms FreeSurfer-formatted ROI names from the
#' DK atlas into labels that are compatible with the `ggseg` R package for plotting
#' data onto a brain surface parcellation
#' 
#' @param .data larger merged dataframe
#' @param feature_roi_col name of existing column containing FreeSurfer-formatted ROI names from the DK atlas
#' @param label_col name of new ggseg label column (defaults to "label")
#' 
#' @return dataframe with new column
#' 
#' @export

add_ggseg_label <- function(.data, feature_roi_col = "feature_roi", label_col = "label") {

    map_ggseg_label <- function(s){
        s_split <- stringr::str_split(s, "\\.")[[1]]
        if (stringr::str_detect(s, "ctx.")) {
            return(paste(s_split[2], s_split[3], sep = "_"))
        } else if (any(stringr::str_detect(s, c("left","right")))) {
            return(paste(stringr::str_to_title(s_split), collapse = "-"))
        } else {
            return(s)
        }
    }

    .data <- .data %>%
        dplyr::mutate("{label_col}" := purrr::map_chr(.data[[feature_roi_col]], map_ggseg_label))

    return(.data)

}

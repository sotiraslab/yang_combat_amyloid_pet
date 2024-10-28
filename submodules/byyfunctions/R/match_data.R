#' @title Match data by date of entry
#' 
#' @description Match a recorded data point to another recorded data point that is closest
#' in time to it
#' 
#' @details For example, match the CDR value closest in time to each subject's
#' PET scan. Each row containing non-NA data in the reference column will be
#' matched with a row containins non-NA data in the matching column which is closest
#' in value to the reference value. User may also specify additional columns to pull
#' data from the matched rows.
#' 
#' @param .data dataframe containing the reference data and matching data
#' @param ref_col name of column to serve as reference value to match to
#' @param match_col name of column to get values to match
#' @param select_col names of column(s) to pull additional data from matched rows
#' @param group_col name of subject identifier column ("RID" for ADNI, "Subject" for
#'     OASIS)
#' @param max_diff (optional) the maximum allowed time difference for a valid match;
#'     if no entires fall within this allowed time difference, NA will be matched
#' @param match_df (optional) a separate dataframe or tibble from which matched data
#'     will be pulled from. If specified, then `match_col` and `select_col` will be
#'     pulled from `match_df`. If not specified, then the aforementioned columns
#'     are pulled from `.data`
#' @param include_idx (optional) if TRUE, then includes an additional column which
#'     has the global index of the matched row for each reference value; note that
#'     this index will be incorrect once the output dataframe is filtered
#' 
#' @return .data with additional columns containing matched data
#' 
#' @section Notes:
#' When using dplyr::group_by, all subsequent functions on a dataframe
#' are performed per group. For our purposes, since we want to find the row index
#' corresponding to the matched entry for each reference value, the returned
#' index is of the group dataframe rather than the entire dataframe, which is not
#' the desired result. To combat this, we make a new column which indicates the
#' global index; that way, our matching function can reference this column and
#' return the global index rather than the group-wise index
#' 
#' @export

# TODO: address error where no matches made ("Logical subscript `match_idx[[1]]` must be size 1 or 507, not 3041.")
# TODO: add argument to control new column names (e.g. whether to have "_match_" prefix)

match_data <- function(
    .data,
    ref_col,
    match_col,
    select_col,
    group_col = "RID",
    max_diff = NULL,
    match_df = NULL,
    include_idx = FALSE
) {

    # ----------------------------
    # ----- helper functions -----
    # ----------------------------

    match_nearest_val_vec <- function(ref_vec, match_list) {

        # Returns the index of vector "match_vec" that is closest in value
        # to each element in the vector "ref_vec"

        match_nearest_val <- function(ref_val) {
            # if ref_val is NA, return NA
            if (is.na(ref_val)) {
                return(NA)
            # else return index of closest value in match_vec (ignore NA in match_vec)
            } else {
                # compute absolute difference
                abs_diff <- abs(ref_val - match_vec)

                # get idx of minimum abs diff
                match_idx <- which.min(abs_diff)

                if (identical(match_idx, integer(0))) {
                    # if match_idx is integer(0) (i.e. no matches), return NA
                    return(NA)
                } else if (!is.null(max_diff) && min(abs_diff, na.rm = TRUE) > max_diff) {
                    # if max_diff specified and no absolute difference falls below this value, return NA
                    return(NA)
                } else {
                    # else return idx
                    return(match_idx)
                }
            }
        }

        match_vec <- match_list[[1]]
        row_idx_global <- match_list[[2]]

        match_idx_vec <- sapply(ref_vec, match_nearest_val)

        # BYY 05/11/23: if match_idx_vec is all NA, then after indexing, the resulting vector
        # will be the same length of row_idx_global, which is *not* the desired behavior; to
        # address this, we check if match_idx_vec is all NA, and if so return an all-NA vector
        # with same length as match_idx_vec
        # - I think this is because when you index with a vector of all NA, R's recycling behavior
        #   kicks in, in which it thinks that you are indexing with only a single NA, thus giving
        #   you a vector which length matches that of the vector being indexed
        # - see https://search.r-project.org/CRAN/refmans/vctrs/html/vector_recycling_rules.html
        if (all(is.na(match_idx_vec))) {
            return(match_idx_vec)
        } else {
            return(row_idx_global[match_idx_vec])
        }

    }
    
    get_group_vector <- function(.data, match_col, group_col, group_id) {

        # helper function to get subset vector from match_df, depending on current group
        # https://stackoverflow.com/questions/28432411/filter-two-data-frames-by-the-same-group-variables-in-dplyr

        group_vector <- .data %>%
            dplyr::filter(.data[[group_col]] %in% unique(group_id)) %>%
            dplyr::select(dplyr::all_of(c(match_col, "row_idx_global")))
        
        group_list <- list(group_vector %>% dplyr::pull(.data[[match_col]]), group_vector$row_idx_global)

        return(group_list)

    }

    # ----------------
    # ----- main -----
    # ----------------

    # get row indices of matched data
    idx_col_name <- paste0(ref_col, "_match_", match_col, "_idx")

    if (is.null(match_df)) {
        match_idx <- .data %>%
            dplyr::select(dplyr::all_of(c(ref_col, match_col, group_col))) %>%
            tibble::rowid_to_column(var = "row_idx_global") %>% # make new column with global idx
            dplyr::group_by(.data[[group_col]]) %>%
            dplyr::mutate("{idx_col_name}" := match_nearest_val_vec(.data[[ref_col]], list(.data[[match_col]], row_idx_global))) %>%
            dplyr::ungroup() %>%
            dplyr::select(dplyr::all_of(idx_col_name))
    } else {
        match_df <- match_df %>% tibble::rowid_to_column(var = "row_idx_global")  # add idx column to match_df
        match_idx <- .data %>%
            dplyr::select(dplyr::all_of(c(ref_col, group_col))) %>%
            dplyr::group_by(.data[[group_col]]) %>%
            dplyr::mutate(
                "{idx_col_name}" := match_nearest_val_vec(
                                        .data[[ref_col]],
                                        match_df %>% get_group_vector(match_col, group_col, .data[[group_col]])
                                    )
            ) %>%
            dplyr::ungroup() %>%
            dplyr::select(dplyr::all_of(idx_col_name))
    }
    
    # get matched columns to select + rename
    select_col <- c(select_col, match_col) # add match col to list of selected columns
    col_rename <- paste0(ref_col, "_match_", select_col) # get renamed column names
    
    # extract data (either from .data or match_df)
    if (is.null(match_df)) {
        match_data <- .data[match_idx[[1]], select_col]
    } else {
        match_data <- match_df[match_idx[[1]], select_col]
    }
    
    # rename columns
    colnames(match_data) <- col_rename 

    # concatenate with main dataframe
    if (include_idx) {
        .data <- dplyr::bind_cols(.data, match_data, match_idx)
    } else {
        .data <- dplyr::bind_cols(.data, match_data)
    }

    return(.data)

}

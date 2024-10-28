#' @title Compute differences between adjacent ordered rows for a dataframe
#' 
#' @description This function first orders the data by the data column, then computes
#' the difference between the current row and the row before it. The differene is
#' stored in a new column called `{col}.diff`. User may also specify grouping variables,
#' such that ordering and difference computation is only done within each group. For the
#' first row in each group, the difference is NA
#'  
#' @param .data dataframe or tibble
#' @param col name of column containing data to order and compute difference
#' @param group_col (optional) name(s) of column(s) for grouping
#' 
#' @return dataframe or tibble with adjacent difference column, named `{col}.diff`
#' 
#' @export
get_diff_adjacent <- function(.data, col = "age", group_col = NULL) {

    df <- .data %>%
        dplyr::group_by(dplyr::pick(dplyr::all_of(group_col))) %>%
        dplyr::arrange(.data[[col]], .by_group = TRUE) %>%
        dplyr::mutate(dplyr::across(
            dplyr::all_of(col),
            ~ .x - lag(.x),
            .names = "{.col}.diff"
        )) %>%
        dplyr::ungroup()

    return(df)

}
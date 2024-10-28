#' @title Rename columns in a dataframe
#' 
#' @description This function renames column given a vector of old column names and
#' corresponding vector of new column names
#
#' @param .data dataframe with columns to rename
#' @param col_old vector of old column names
#' @param col_new vector of new column names; must be the same
#'     length as col_old
#
#' @return dataframe with renamed columns
#' 
#' @export
rename_columns <-
function(.data, col_old, col_new) {

    # Rename column given a vector of old column names and
    # corresponding new column names
    #
    # Parameters
    # ----------
    #   .data :     dataframe with columns to rename
    #   col_old :   vector of old column names
    #   col_new :   vector of new column names; must be the same
    #             length as col_old
    #
    # Returns
    # -------
    #   dataframe with renamed columns

    .data <- .data %>%
        dplyr::rename_with(~ col_new[which(col_old == .x)], .cols = col_old)

    return(.data)

}

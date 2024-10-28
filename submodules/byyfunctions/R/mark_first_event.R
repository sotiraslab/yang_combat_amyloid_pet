#' @title Mark where a partcular event occured first within each group
#' 
#' @description This function creates a new boolean column (named with argument
#' `new_col`) which indicates where a particular event occured first within a
#' single group (usually a subject)
#' 
#' @details For this function to operate, one column must contain date or time of
#' event information (ex. date of a particular AV45 PET scan), and another column
#' must contain unique identifiers of groups (ex. subject IDs). These two are
#' inputted into the function as arguments `date_col` and `group_col`, respectively.
#' For each group ID, there will be only one TRUE entry which corresponds to the
#' earliest date in column `date_col`
#' 
#' @param .data larger merged dataframe
#' @param date_col column containing event dates
#' @param group_col column containing group identifiers (usually subject IDs)
#' @param new_col name of new column which indicates first events
#' 
#' @return dataframe with new column
#' 
#' @export
mark_first_event <-
function(.data, date_col, group_col, new_col) {
    .data <- .data %>%
        dplyr::group_by(.data[[group_col]]) %>%
        dplyr::mutate(
            "{new_col}" := if (all(is.na(.data[[date_col]]))) {NA} else {.data[[date_col]] == min(.data[[date_col]], na.rm = TRUE)}
        ) %>%
        dplyr::ungroup(.data[[group_col]])

    return(.data)
}

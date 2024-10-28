#' @title Add suffix to column names
#' 
#' @description This function appends a suffix to all columns of a dataframe or tibble.
#' One may additionally exclude columns from being appended by specifying them in the
#' `exclude` argument
#' 
#' @param .data dataframe or tibble
#' @param suffix suffix to append to column names
#' @param exclude character vector containing columns to exclude; if not specified, then
#'     all columns will be appended by default
#' 
#' @return dataframe or tibble with new column names
#' 
#' @section Notes:
#' - see https://stackoverflow.com/questions/52966616/exclude-some-columns-from-adding-a-suffix-in-r
#'   for more details
#' 
#' @export
colname_add_suffix <-
function(.data, suffix, exclude=NULL) {

    if (!is.null(exclude)) {
        # add suffix to all columns in .data, EXCLUDING columns listed in exclude
        # https://stackoverflow.com/questions/52966616/exclude-some-columns-from-adding-a-suffix-in-r
        col_idx <- which(colnames(.data) %in% exclude)
        colnames(.data)[-col_idx] <- paste(colnames(.data)[-col_idx], suffix, sep = "")
    }
    else {
        colnames(.data) <- paste(colnames(.data), suffix, sep = "")
    }
    
    return(.data)

}

#' @title Load multiple CSVs in a single directory
#' 
#' @description This function loads multiple CSVs from a specified directory and
#' stores them in either a list (bind = FALSE) or a single dataframe (bind = TRUE)
#'  
#' @param d name of directory containing CSVs
#' @param vector if TRUE, then load CSVs without header; this is the case when the
#' CSVs only store a single vector of data
#' @param bind if TRUE, then load CSVs and bind into one dataframe, otherwise store
#' each individual dataframe in a list
#' 
#' @return either list of dataframes or a single dataframe
#' 
#' @export
load_csv_dir <- function(d, vector = FALSE, bind = FALSE) {

    p <- list.files(
        path = d,
        pattern = "*.csv",
        full.names = TRUE
    )

    if (bind) {

        l <- readr::read_csv(p, show_col_types = FALSE, id = "filename")

    } else {
    
        if (vector) {
            f <- ~ readr::read_csv(.x, show_col_types = FALSE, col_names = FALSE) %>% dplyr::pull(1)
        } else {
            f <- ~ readr::read_csv(.x, show_col_types = FALSE)
        }

        l <- purrr::map(p, f)
        names(l) <- p %>% basename %>% stringr::str_remove(".csv")
        
    }

    return(l)

}
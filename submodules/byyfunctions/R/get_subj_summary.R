#' @title Get subject-specific dataframe from a larger merged dataframe
#' 
#' @description Return dataframe where each row is a single subject, and each columns is a
#' subject-level field (e.g. date of first AV45+ scan, or first CDR record)
#' 
#' @details The purpose of this function is to return a dataframe where each row is on
#' the subject level, rather than on the entry level (as is the case with the
#' main ADNIMERGE tables). Since all rows will have the same value within-subject, the first row
#' of each subject is used to populate the subject-specific dataframe
#' 
#' @param .data larger merged dataframe
#' @param column_names vector of columns names to extract
#' @param group_name name of column to group entries by (usually subject ID)
#' 
#' @return dataframe on the subject level with specified columns selected
#' 
#' @export
get_subj_summary <-
function(.data, column_names, group_name="Subject") {
  
    # Return dataframe where each row is a single subject, and each columns is a
    # subject-level field (e.g. date of first AV45+ scan, or first CDR record);
    # the purpose of this function is to return a dataframe where each row is on
    # the subject level, rather than on the entry level (as is the case with the
    # main ADNIMERGE tables)
    # 
    # NOTE: since all rows will have the same value within-subject, the first row
    # of each subject is chosen
    # 
    # Parameters
    # ----------
    #   .data :        ADNI merged dataframe
    #   column_names : vector of columns names to extract
    #   group_name :   name of column to group entries by (usually subject ID)
    # 
    # Returns
    # -------
    #   dataframe on the subject level with specified columns selected
  
    .data <- .data %>%
        dplyr::group_by(.data[[group_name]]) %>% # group by Subject
        dplyr::slice(1) %>% # select the *first* row
        dplyr::select(dplyr::all_of(c(column_names, group_name))) # select columns from string vector

    return(.data)
  
}

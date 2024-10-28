#' @title Concatenate strings while propagating NAs
#' 
#' @description This function acts exactly like "paste0", only it propagates
#' NAs. This is useful for vectors which contain NAs, but the desired output
#' is NA rather than "NA{suffix}"
#' 
#' @param s_vec character vector
#' @param s_paste string to paste
#' 
#' @return pasted character vector with NAs propagated
#' 
#' @export
paste0_na <-
function(s_vec, s_paste) {
    return(ifelse(is.na(s_vec), NA, paste0(s_vec, s_paste)))
}

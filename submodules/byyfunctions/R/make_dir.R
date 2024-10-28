#' @title Recursively make directory if doesn't exist
#' 
#' @description This function first checks to see whether a directory exists. If 
#' it doesn't, it creates the directory and any necessary subdirectories. If it
#' does exist, it does nothing
#' 
#' @param d path to directory
#' @param verbose if TRUE, then print status to terminal
#' 
#' @export
make_dir <-
function(d, verbose = FALSE) {
    if (!dir.exists(d)){
        if (verbose) {print(paste("++ creating", d, "++", sep = " "))}
        dir.create(d, recursive = TRUE)
    } else {
        if (verbose) {print(paste("++", d, "already exists ++", sep = " "))}
    }
}

#' @title Column names of FreeSurfer ROIs in ADNI and OASIS
#' 
#' @description A tidy tibble containing column labels of FreeSurfer
#' ROI-based features (both amyloid SUVR and volumetric) from the ADNI
#' and OASIS datasets, along with corresponding FreeSurfer labels, as
#' denoted in the color look-up table. Note that column names also
#' contain suffixes which roughly indicate which original spreadsheets
#' the columns came from (for example ".pup" in OASIS means the column
#' came from the PUP spreadsheet)
#' 
#' @docType data
#' 
#' @usage data(roi_features)
#' 
#' @format a tibble
#' 
"roi_features"
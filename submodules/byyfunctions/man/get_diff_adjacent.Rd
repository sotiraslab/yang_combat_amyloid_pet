% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_diff_adjacent.R
\name{get_diff_adjacent}
\alias{get_diff_adjacent}
\title{Compute differences between adjacent ordered rows for a dataframe}
\usage{
get_diff_adjacent(.data, col = "age", group_col = NULL)
}
\arguments{
\item{.data}{dataframe or tibble}

\item{col}{name of column containing data to order and compute difference}

\item{group_col}{(optional) name(s) of column(s) for grouping}
}
\value{
dataframe or tibble with adjacent difference column, named \verb{\{col\}.diff}
}
\description{
This function first orders the data by the data column, then computes
the difference between the current row and the row before it. The differene is
stored in a new column called \verb{\{col\}.diff}. User may also specify grouping variables,
such that ordering and difference computation is only done within each group. For the
first row in each group, the difference is NA
}

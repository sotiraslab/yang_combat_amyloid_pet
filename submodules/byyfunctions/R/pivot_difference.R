#' @title Compute difference in long-formatted values using pivot_wider
#' 
#' @description This function uses a series of tidyr pivotting operations to compute
#' the difference in values compared to a reference name. This requires that the df
#' containing the data is in long format, and that the names of the variable types
#' are stored in a separate column.
#' 
#' @details First, the function runs `tidyr::pivot_wider` to convert the dataframe
#' into wide format. Then it computes the difference between each widened column and
#' the reference column. Finally it runs `tidyr::pivot_longer` to convert the dataframe
#' back into long format. Differences in values are stored in the original values column,
#' which by default has its values replaced by the differences.
#' 
#' @param .data dataframe with data to compute difference
#' @param names_col name of col containing names of variable groupings
#' @param values_col name(s) of col containing values; can accept multiple value columns
#' @param ref name of the variable group to serve as the reference; must exist within the
#' `values_col` column
#' @param id_col (optional) vector of column names to serve as ID cols for pivot_wider
#' 
#' @return dataframe with differences in values
#' 
#' @export
pivot_difference <- function(.data, names_col, values_col, ref, id_col = NULL) {

    # TODO: look into implementing tidyselect interface for `id_col`
    #   - https://tidyselect.r-lib.org/articles/tidyselect.html

    names_list <- .data %>%
        dplyr::pull(names_col) %>%
        unique

    # first pivot longer in case of multiple values
    df <- .data %>%
        tidyr::pivot_longer(
            cols = starts_with(values_col),
            names_to = "val_name",
            values_to = "val_val"
        )

    # pivot wider
    if (is.null(id_col)) {
        df <- df %>%
            tidyr::pivot_wider(
                names_from = dplyr::all_of(names_col),
                values_from = val_val
            )
    } else {
        df <- df %>%
            tidyr::pivot_wider(
                id_cols = dplyr::all_of(id_col),
                names_from = dplyr::all_of(names_col),
                values_from = val_val
            )
    }

    # compute difference, pivot longer, then pivot wider
    df <- df %>%
        dplyr::mutate(dplyr::across(
            dplyr::all_of(names_list),
            ~ .x - .data[[ref]],
            .names = "diff__{.col}"
        )) %>%
        tidyr::pivot_longer(
            cols = dplyr::contains(names_list),
            names_pattern = "(diff__|)(.*)",
            names_to = c("diff", names_col),
            values_to = "val_val"
        ) %>%
        dplyr::mutate(
            val_name = dplyr::if_else(
                diff == "diff__",
                stringr::str_c(val_name, ".diff"),
                val_name
            )
        ) %>%
        dplyr::select(-diff) %>%
        tidyr::pivot_wider(
            names_from = val_name,
            values_from = val_val
        )
    
    return(df)

}

# alt
# pivot_wider(
#     id_cols = -magnitude,
#     names_from = harmonization_method,
#     values_from = effsize
# ) %>%
# mutate(across(
#     c(cl, starts_with("combat")),
#     ~ .x - suvr,
#     .names = "diff.{.col}"
# )) %>%
# pivot_longer(
#     cols = c(cl, starts_with("combat"), starts_with("diff")),
#     names_pattern = "(diff\\.|)(.*)",
#     names_to = c("diff", "harmonization_method"),
#     values_to = c("effsize")
# ) %>%
# mutate(diff = if_else(diff == "diff.", "diff", "val"))
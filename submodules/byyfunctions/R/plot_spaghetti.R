#' @title Create spaghetti plot for visualizing longitudinal subject-level data
#' 
#' @description This function plots a spaghetti plot to visualize longitudinal data,
#' where the x-axis is the time axis and the y-axis is a discrete subject index
#' 
#' @details This function requires a dataframe in long-ish format where each row
#' corresponds to a record at a single time point for a single subject. The table
#' at minimum must contain a subject ID column and a time column. Values of other
#' features (e.g. intracranial volume, PET SUVR, amyloid-positivity, etc.) are
#' encoded via the shape (for discrete values) and color (for continuous values)
#' aesthetics. Additionally, one may order the y-axis indices based on `col_order`,
#' which contains subject-level values (i.e. a single value corresponding to each
#' subject)
#' 
#' @param .data dataframe or tibble containing longitudinal data
#' @param col_x column containing time values; this will typically be a
#' longitudinal metric, such as age or days from enrollment
#' @param col_group column containing unique subject identifiers
#' @param col_color (optional) column containing values to encode point colors
#' @param col_shape (optional) column containins values to encode point shapes
#' @param col_order (optional) column used to determine the ordering of the
#' y-axis index. For example, one may choose to order each "spaghetti string"
#' by date of AV45-positivity rather than age at enrollment. If not specified,
#' y-axis indices will be ordered by the minimum `col_x` value by default
#' 
#' @return ggplot spaghetti plot
#' 
#' @export
plot_spaghetti <- function(.data, col_x, col_group, col_color=NULL, col_shape=NULL, col_order=NULL, dodge=NULL) {

    # select relevant columns
    col_to_select <- c(col_group, col_x)
    for (col in c(col_color, col_shape, col_order)) {
        col_to_select <- if (!is.null(col)) c(col_to_select, col) else col_to_select
    }
    new_df <- .data %>%
        dplyr::select(all_of(col_to_select)) %>% # select columns from string vector
        tidyr::drop_na() # filter out NA

    # sort by minimum col_order or minimum col_x, assign order index
    if (!is.null(col_order)) {
        sort_idx_df <- new_df %>%
            dplyr::group_by(.data[[col_group]]) %>% # group by RID
            dplyr::summarise(min_x = min(.data[[col_order]])) # compute each group's minimum col_order value
    } else {
        sort_idx_df <- new_df %>%
            dplyr::group_by(.data[[col_group]]) %>% # group by RID
            dplyr::summarise(min_x = min(.data[[col_x]])) # compute each group's minimum col_x value
    }

    # sort minimum values and assign order index
    sort_idx_df <- sort_idx_df[order(sort_idx_df$min_x),] # sort x values
    sort_idx_df$idx <- as.numeric(row.names(sort_idx_df)) # assign order index

    # merge with new_df
    new_df <- dplyr::left_join(new_df, sort_idx_df, by=col_group) # merge with new_df

    # define geom_point
    my_geom_point <-ggplot2::geom_point(
        size=2,
        ggplot2::aes(
            x=.data[[col_x]],
            group=if(is.null(dodge)) NULL else if(dodge == "color") as.factor(.data[[col_color]]) else if(dodge == "shape") .data[[col_shape]] else NULL,
            color=if(is.null(col_color)) NULL else as.factor(.data[[col_color]]),
            shape=if(is.null(col_shape)) NULL else .data[[col_shape]]
        ),
        position = if(is.null(dodge)) "identity" else ggstance::position_dodgev(height = 0.5)
    )

    # plot spaghetti plot
    spaghetti_plot <- new_df %>% ggplot2::ggplot(aes(y = idx)) +
        ggplot2::geom_line(mapping = ggplot2::aes(x = .data[[col_x]], group = .data[[col_group]]), linewidth=1) +
        my_geom_point +
        ggplot2::geom_line(color = "gray40", mapping = ggplot2::aes(x = min_x)) +
        ggplot2::labs(
            x = col_x,
            y = "Subject",
            color = col_color,
            shape = col_shape
        )
  
    return(spaghetti_plot)
  
}

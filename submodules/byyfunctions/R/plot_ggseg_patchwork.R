#' @title Plot multiple ggseg plots
#' 
#' @description This function plots multiple brain surface + aseg ggseg plots, one per
#' group, then stores them in a list. This is useful for combining multiple plots into
#' one plot via patchwork's `wrap_plots` function.
#' 
#' @param .data dataframe or tibble containing data to plot on surface + aseg
#' @param fill_col string of column containing fill data
#' @param group_col string of column containing the group labels; this **must** be a factor
#' column
#' @param scale ggplot2 fill scale; this will be applied to *every* plot generated
#' @param cb_title (optional) title of colorbar; if none specified, defaults to `fill_col`
#' @param aseg_side string to indicate the view to plot for aseg if you only 
#' want to plot one; either 'coronal' or 'sagittal'
#' 
#' @return list of ggseg plots, with names in the form of "{group_level}.{surf|aseg}"
#' 
#' @export
plot_ggseg_patchwork <- function(.data, fill_col, group_col, scale, cb_title = NULL, aseg_side = NULL) {

    if (is.null(cb_title)) {
        cb_title = fill_col
    }

    cb <- ggplot2::guides(
        fill = ggplot2::guide_colorbar(
            barheight = ggplot2::unit(2.5, "inch"),
            title = cb_title
        )
    )

    ggseg_list <- list()
    for (g in .data[[group_col]] %>% levels) {

        df <- .data %>% dplyr::filter(.data[[group_col]] == g)
            
        ggseg_list[[paste0(g, ".surf")]] <- df %>%
            byyfunctions::plot_ggseg(fill_col = fill_col, aseg_bool = FALSE, aseg_side = aseg_side) +
            ggplot2::labs(title = g) +
            ggplot2::theme(plot.title = ggplot2::element_text(vjust = 0, hjust = 0)) +
            scale + cb
            
        ggseg_list[[paste0(g, ".aseg")]] <- df %>%
            byyfunctions::plot_ggseg(fill_col = fill_col, aseg_bool = TRUE, aseg_side = aseg_side) +
            scale + cb
    }
    
    return(ggseg_list)

}
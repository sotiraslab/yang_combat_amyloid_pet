#' @title Plot ggseg surface + aseg plot
#' 
#' @description This function plots both a brain surface plot and brain aseg (volumetric)
#' plot using ggseg, then combines them using patchwork. The colorbars are collected
#' between the two plots, so the user should only specify one fill scale for the output
#' patchwork plot
#' 
#' @param .data dataframe or tibble containing data to plot on surface + aseg
#' @param fill_col string of column containing fill data
#' @param title_str (optional) title of plot
#' @param cb_title (optional) title of colorbar; if none specified, defaults to `fill_col`
#' @param aseg_side string to indicate the view to plot for aseg if you only 
#' want to plot one; either 'coronal' or 'sagittal'
#' 
#' @return patchwork plot object
#' 
#' @export
plot_surf_aseg <- function(.data, fill_col, title_str = NULL, cb_title = NULL, aseg_side = NULL) {

    if (is.null(cb_title)) {
        cb_title = fill_col
    }
    
    cb <- ggplot2::guides(fill = ggplot2::guide_colorbar(title = cb_title))

    surf <- .data %>%
        byyfunctions::plot_ggseg(fill_col = fill_col, aseg_bool = FALSE, aseg_side = aseg_side) +
        cb
    aseg <- .data %>%
        byyfunctions::plot_ggseg(fill_col = fill_col, aseg_bool = TRUE, aseg_side = aseg_side) +
        cb
    ggseg_layout <- "
    A#
    AB
    A#
    "
    surf_aseg <- (surf | aseg) +
        patchwork::plot_layout(design = ggseg_layout, guides = "collect") &
        patchwork::plot_annotation(
            title = title_str,
            theme = ggplot2::theme(
                plot.title = ggplot2::element_text(vjust = -100, hjust = 0.5)
            )
        )
    
    return(surf_aseg)

}
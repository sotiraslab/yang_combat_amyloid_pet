#' @title Plot data on a brain surface or coronal slice using ggseg
#' 
#' @description This function plots continuous data on a brain surface or coronal
#' slice using the package `ggseg`
#' 
#' @details This function assumes that the input dataframe has a column named
#' "label" which contains ggseg-specific ROI labels of the FreeSurfer DK atlas
#' 
#' @param .data dataframe or tibble containing data to plot
#' @param fill_col name of column containing continuous-valued data
#' @param aseg_bool if true, plot on coronal slice instead of brain surface;
#' useful for visualizing subcortical ROIs
#' @param aseg_side string to indicate the view to plot for aseg if you only 
#' want to plot one; either 'coronal' or 'sagittal'
#' 
#' @return ggseg plot
#' 
#' @export
plot_ggseg <- function(.data, fill_col, aseg_bool = FALSE, aseg_side = NULL) {

    # TODO: add argument for faceting
    
    ggseg_plot <- .data %>% ggplot2::ggplot()
    
    if (aseg_bool) {
        ggseg_plot <- ggseg_plot + ggseg::geom_brain(
            atlas = ggseg::aseg,
            mapping = ggplot2::aes(fill = .data[[fill_col]]),
            side = aseg_side
        )
    } else {
        ggseg_plot <- ggseg_plot + ggseg::geom_brain(
            atlas = ggseg::dk,
            position = ggseg::position_brain(hemi ~ side),
            mapping = ggplot2::aes(fill = .data[[fill_col]])
        )
    }
    
    ggseg_plot <- ggseg_plot +
        ggplot2::labs(fill = fill_col) +
        ggplot2::theme_void()
    
    return(ggseg_plot)
    
}
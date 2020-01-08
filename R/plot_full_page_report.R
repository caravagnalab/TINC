#' Plot a full page detailed report of the TINC analysis
#'
#' @description The report is a multi-panel figure that contains the following
#' information. TODO fill in.
#'
#' @param x A TINC analysis computed with \code{autofit}.
#'
#' @return A multi-panle \code{ggplot} figure.
#' @export
#'
#' @examples
#' plot_full_page_report(autofit(random_TIN(), FAST = TRUE))
plot_full_page_report = function(x)
{

  if(!inherits(x, "tin_obj")) stop("Not a TINC object .... run autofit(.) first, aborting.")

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Report assembly
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  # Plot input data
  dataset_plot_raw = TINC:::plot_raw(dataset = x)

  # Plots of all the data
  data_fit_panel = cowplot::plot_grid(
    dataset_plot_raw,
    TINC:::plot_sample_contamination(x),
    rel_widths = c(1, .7),
    align = 'h',
    labels = c("A", "B"),
    nrow = 1,
    ncol = 2
  )

  # data_fit_panel =
  #   ggpubr::annotate_figure(data_fit_panel,
  #                           top = ggpubr::text_grob(
  #                             paste0('Input: n = ', nrow(file, '\n'),
  #                             hjust = 0,
  #                             x = 0,
  #                             size = 18
  #                           ))

  fit_panel =
    cowplot::plot_grid(
      x$fit$mobster_analysis$plot,
      x$fit$BMix_analysis$plot,
      nrow = 1,
      ncol = 2,
      align = 'h',
      axis = 'bt',
      labels = c('C', 'D')
    )

  sq_panel = cowplot::plot_grid(
    TINC:::plot_contamination_full_size(x),
    TINC:::plot_contamination_zoom(x),
    x$fit$VIBER_analysis$plot,
    nrow = 1,
    ncol = 3,
    align = 'h',
    axis = 'bt',
    rel_widths = c(2, 1, 1),
    labels = c('E', '', 'D')
  )

  cowplot::plot_grid(
    data_fit_panel,
    fit_panel,
    sq_panel,
    nrow = 3,
    ncol = 1,
    align = 'v',
    axis = 'lr',
    rel_heights = c(.8, .7, .7)
  )
}

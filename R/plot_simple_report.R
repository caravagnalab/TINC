#' Plot a strip short report of the TINC analysis.
#'
#' @description The report is a multi-panel figure that contains the following
#' information. TODO fill in.
#'
#' @param x A TINC analysis computed with \code{autofit}.
#'
#' @return A multi-panel \code{ggplot} figure.
#' @export
#' @importFrom ggplot2 geom_point guides scale_fill_manual coord_cartesian scale_fill_brewer 
#' geom_raster geom_abline geom_label coord_polar facet_wrap geom_bar coord_flip
#' scale_x_continuous scale_fill_viridis_c guide_colorbar element_blank guide_legend
#' element_rect alpha scale_color_manual
#' @examples 
#' data(fit_example, package = 'TINC')
#' plot_simple_report(fit_example)
plot_simple_report = function(x)
{
  if(!inherits(x, "tin_obj")) stop("Not a TINC object .... run autofit(.) first, aborting.")

  # Labels
  mut_load = sum(x$data$OK_tumour)
  all_clonal = mobster::Clusters(x$fit$mobster_analysis$fit$best) %>%
    dplyr::filter(cluster == x$fit$mobster_analysis$clonal_cluster) %>%
    nrow

  used_clonal = x$fit$mobster_analysis$clonal_mutations %>% length
  p_clonal = round((used_clonal/all_clonal) * 100)

  label = paste(
    paste0("Tumour mutational burden: n = ", mut_load, " (total)"), '\n',
    paste0("        Clonal mutations: n = ", all_clonal, " (", used_clonal, ' ~ ', p_clonal, "% with confidence >", x$params$cutoff_lv_assignment, ")")
  )

  # data_plot = plot_raw(x$fit)
  stats_plot = plot_sample_contamination(x, assemble = F)
  cont_plot = plot_contamination_full_size(x) +
    labs(
      title = bquote(bold("Summary")),
      subtitle =  paste0("Clonal mutations: n = ",
                         all_clonal, "\n(n = ", used_clonal, ' ~ ',
                         p_clonal, "% with confidence >", x$params$cutoff_lv_assignment, ")"),
      caption = paste0("Tumour mutational burden: n = ", mut_load, " (total)\nFAST analysis: ", x$params$FAST)
    )

  # Sample pie
  # sample_pies = ggpubr::ggarrange(
  #   stats_plot$normal,
  #   stats_plot$tumour,
  #   ncol = 2,
  #   nrow = 1,
  #   common.legend = TRUE,
  #   legend = 'bottom')

  figure = cowplot::plot_grid(
    cont_plot,
    stats_plot,
    ncol = 2,
    nrow = 1,
    align = 'h',
    axis = 'bt'
  )

  return(figure)

#
#   figure =
#     ggpubr::annotate_figure(
#       figure,
#       bottom =
#         ggpubr::text_grob(label, just = 'left', x = 0)
#     )


}

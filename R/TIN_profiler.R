# Assemble results
TIN_profiler = function(dataset)
{
  x %>%
    full_join(
      dataset$mobster_analysis$output %>% dplyr::select(id, tumour.cluster),
      by = c("id")
    ) %>%
    full_join(
      dataset$BMix_analysis$output %>% dplyr::select(id, normal.cluster),
      by = c("id")
    )

  dataset$joint = dataset$joint %>%
    full_join(
      dataset$mobster_analysis$output %>% dplyr::select(id, tumour.cluster),
      by = c("id")
      ) %>%
    full_join(
      dataset$BMix_analysis$output %>% dplyr::select(id, normal.cluster),
      by = c("id")
    )

  is_cl_tum = dataset$joint$tumour.cluster %in% dataset$mobster_analysis$clonal_cluster
  is_cl_nor = dataset$joint$normal.cluster %in% dataset$BMix_analysis$clonal_cluster

  dataset$joint$is_clonal = is_cl_tum & is_cl_nor

  dataset$joint = dataset$joint %>%
    dplyr::select(
      chr,
      from,
      to,
      ref,
      alt,
      id,
      ends_with('count'),
      ends_with('tumour'),
      ends_with('normal'),
      ends_with('cluster'),
      used,
      is_clonal
    )

  dataset
}

  # TIT = dataset$mobster_analysis$estimated_purity / 2
  # TIN = dataset$$estimated_purity / 2

  # ggplot(
  #   data_table %>% mutate(is_clonal = ifelse(is_clonal, 'Clonal', 'Subclonal')),
  #        aes(x = VAF.germline, y = VAF.tumour, color = is_clonal)) +
  #   geom_point(
  #     size = .3,
  #     alpha = 0.1,
  #     color = 'black'
  #   ) +
  #   # scale_x_continuous(limits=c(0, 1)) +
  #   facet_zoom(xlim = c(0, 0.05)) +
  #   ylim(0, 1) +
  #   mobster:::my_ggplot_theme()

#
#   figure = ggplot(data_table %>% mutate(is_clonal = ifelse(is_clonal, 'Clonal', 'Subclonal or CNA')),
#                   aes(x = VAF.normal, y = VAF.tumour, color = is_clonal)) +
#     geom_point(
#       data = dataset$joint,
#       size = .3,
#       alpha = 0.1,
#       color = 'black'
#     ) +
#     geom_point(size = .7, alpha = .3) +
#     geom_abline(linetype = 'dashed',
#                 color = 'red',
#                 size = .3) +
#     geom_hline(
#       yintercept = TIT,
#       linetype = 'dashed',
#       color = 'black',
#       size = .3
#     ) +
#     geom_vline(
#       xintercept = TIN,
#       linetype = 'dashed',
#       color = 'black',
#       size = .3
#     ) +
#     geom_point(
#       data = data.frame(x = TIN, y = TIT),
#       aes(x, y),
#       color = 'black',
#       fill = 'white',
#       size = 3,
#       shape = 21
#     ) +
#     xlim(0, 1) +
#     ylim(0, 1) +
#     guides(color = guide_legend('')) +
#     mobster:::my_ggplot_theme() +
#     labs(
#       title = paste0("Contamination"),
#       x = 'VAF germline',
#       y = "VAF tumour"
#       # subtitle = paste0(
#       #   "Cancer cells in normal sample: ",
#       #   round(TIN * 200, 3),
#       #   '%\n',
#       #   "Cancer cells in tumour sample: ",
#       #   round(TIT * 200, 3),
#       #   '%'
#       # )
#       # caption = paste0('S = ', S, '(', percentage * 100, '%)')
#     ) +
#     scale_color_manual(values = c(`Clonal` = 'cyan4', `Subclonal or CNA` = 'black'))
#
#   # geom_mark_ellipse(aes(fill = is_clonal, label = is_clonal))
#
#
#   data_table = data_table %>%
#     mutate(
#       group = paste(tumour.cluster, '~', germline.cluster)
#     )
#
#   group_labels = data_table$group %>% unique
#   group_labels_colors = pio:::nmfy(group_labels, rep("black", length(group_labels)))
#
#   for(cl in MOBSTER_analysis_results$clonal_cluster)
#     group_labels_colors[grepl(group_labels, pattern = cl)] = 'cyan4'
#
#
#   f2 = ggplot(data_table,
#               aes(
#                 x = VAF.germline,
#                 y = VAF.tumour,
#                 shape = group,
#                 color = group
#               )) +
#     # geom_point(
#     #   data = joint,
#     #   aes(x = VAF.germline,
#     #       y = VAF.tumour),
#     #   size = .3,
#     #   alpha = 0.1,
#     #   color = 'black',
#     #   inherit.aes = FALSE
#     #
#     # ) +
#     # geom_density_2d(data = data_table %>% filter(VAF.germline > 0), size = .1) +
#   geom_point(size = 2, alpha = .7) +
#     geom_abline(linetype = 'dashed',
#                 color = 'red',
#                 size = .3) +
#     geom_hline(
#       yintercept = TIT,
#       linetype = 'dashed',
#       color = 'black',
#       size = .3
#     ) +
#     geom_vline(
#       xintercept = TIN,
#       linetype = 'dashed',
#       color = 'black',
#       size = .3
#     ) +
#     geom_point(
#       data = data.frame(x = TIN, y = TIT),
#       aes(x, y),
#       color = 'black',
#       fill = 'white',
#       size = 3,
#       shape = 21
#     ) +
#     guides(shape = guide_legend("Clusters"), color = FALSE) +
#     mobster:::my_ggplot_theme() +
#     labs(
#       title = paste0(""),
#       x = 'VAF germline',
#       y = "VAF tumour"
#       # subtitle = paste0(
#       #   "Clonal cluster (tumour) : ",
#       #   paste(MOBSTER_analysis_results$clonal_cluster, collapse = ', '),
#       #   " (n = ", length(MOBSTER_analysis_results$clonal_mutations), ')',
#       #   '\n',
#       #   "Germline weigth : ",
#       #   paste(germline_analysis_results$clonal_cluster, collapse = ', ')
#       # )
#     ) +
#     scale_color_manual(values = group_labels_colors)
#
#   vtin = data.frame(
#     variable = c('Tumour cells', 'Normal cells', 'Tumour cells', 'Normal cells'),
#     value = c(TIN * 2, 1 - (TIN * 2), TIT * 2, 1 - (TIT * 2)),
#     sample = c(
#       'Normal sample',
#       'Normal sample',
#       'Tumour sample',
#       'Tumour sample'
#     )
#   )
#
#   f3 = ggplot(data = vtin,
#               aes(x='1', y = value, fill = variable)) +
#     geom_bar(stat = 'identity') +
#     mobster:::my_ggplot_theme() +
#     facet_wrap(~sample, ncol = 2, nrow = 1) +
#     # coord_flip() +
#     coord_polar("y", start = 0, clip = 'off') +
#     labs(
#       title = paste0("Sample composition"),
#       x = element_blank(),
#       y = "Percentage",
#       subtitle = paste0(
#         "TIN: ",
#         round(vtin$value[1] * 100, 2),
#         '%',
#         " (",
#         round(vtin$value[2] * 100, 2),
#         '% germline)\n',
#         "TIT : ",
#         round(vtin$value[3] * 100, 2),
#         '% (',
#         round(vtin$value[4] * 100, 2),
#         '% germline)'
#       )
#     ) +
#     scale_fill_manual(values = c(`Tumour cells` = 'plum4', `Normal cells` = 'plum2')) +
#     guides(fill = guide_legend(''))
#
#
#   return(list(data_table = data_table,
#               summary_plot = figure,
#               zoom_plot = f2,
#               sample_composition = f3,
#               TIN = TIN, TIT = TIT))

  #
  # figure = cowplot::plot_grid(
  #   # plot_input,
  #   figure,
  #   f2,
  #   f3,
  #   ncol = 3,
  #   nrow = 1,
  #   align = 'h', axis = 'x'
  # )
  #
  #
  # return(list(data_table = data_table, plot = figure, TIN = TIN, TIT = TIT))
  #



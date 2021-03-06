plot_contamination_full_size = function(dataset)
{
  TIT = dataset$TIT / 2
  TIN = dataset$TIN / 2
  
  ggplot(
    as_joint(dataset$data) %>%
      mutate(
        is_clonal = ifelse(OK_clonal, 'Clonal', 'Other'),
        is_clonal = ifelse(t_VAF > 0, is_clonal, 'Excluded')
      ),
    aes(x = n_VAF, y = t_VAF, color = is_clonal)
  ) +
    geom_point(size = .3,
               alpha = 0.1,
               color = 'black') +
    geom_point(size = .7, alpha = .3) +
    geom_abline(linetype = 'dashed',
                color = 'red',
                size = .3) +
    geom_hline(
      yintercept = TIT,
      linetype = 'dashed',
      color = 'black',
      size = .3
    ) +
    geom_vline(
      xintercept = TIN,
      linetype = 'dashed',
      color = 'black',
      size = .3
    ) +
    geom_point(
      data = data.frame(x = TIN, y = TIT),
      aes(x, y),
      color = 'black',
      fill = 'white',
      size = 3,
      shape = 21
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    guides(color = guide_legend('', override.aes = list(alpha = 1, size = 2))) +
    my_ggplot_theme() +
    labs(
      title = paste0("Contamination of clonal mutations"),
      x = 'VAF germline',
      y = "VAF tumour"
    ) +
    scale_color_manual(values = c(
      `Clonal` = 'cyan4',
      `Other` = 'black',
      `Excluded` = 'gray'
    )) +
    theme(legend.position = c(0.8, 0.2))
}

plot_contamination_zoom = function(dataset)
{
  TIT = dataset$TIT / 2
  TIN = dataset$TIN / 2
  
  data_table = as_joint(dataset$data) %>%
    mutate(group = paste(tumour.cluster, '~', normal.cluster))
  
  group_labels = data_table$group %>% unique
  group_labels_colors = pio:::nmfy(group_labels, rep("black", length(group_labels)))
  
  for (cl in dataset$fit$mobster_analysis$clonal_cluster)
    group_labels_colors[grepl(group_labels, pattern = cl)] = 'cyan4'
  
  
  ggplot(data_table,
         aes(
           x = n_VAF,
           y = t_VAF,
           shape = group,
           color = group
         )) +
    geom_point(size = 2, alpha = .7) +
    geom_abline(linetype = 'dashed',
                color = 'red',
                size = .3) +
    geom_hline(
      yintercept = TIT,
      linetype = 'dashed',
      color = 'black',
      size = .3
    ) +
    geom_vline(
      xintercept = TIN,
      linetype = 'dashed',
      color = 'black',
      size = .3
    ) +
    geom_point(
      data = data.frame(x = TIN, y = TIT),
      aes(x, y),
      color = 'black',
      fill = 'white',
      size = 3,
      shape = 21
    ) +
    guides(shape = guide_legend("Clusters"), color = FALSE) +
    my_ggplot_theme() +
    labs(title = paste0(""),
         x = 'VAF germline (zoom)',
         y = "VAF tumour (zoom)") +
    scale_color_manual(values = group_labels_colors)
}

plot_sample_contamination = function(dataset, assemble = TRUE)
{
  # fot the plots we half the estimates
  TIT = dataset$TIT
  TIN = dataset$TIN
  
  vtin = data.frame(
    variable = c('Tumour cells', 'Normal cells', 'Tumour cells', 'Normal cells'),
    value = c(TIN, 1 - (TIN), TIT, 1 - (TIT)),
    sample = c(
      'Normal sample',
      'Normal sample',
      'Tumour sample',
      'Tumour sample'
    )
  ) %>%
    group_by(sample) %>%
    mutate(lab.ypos = cumsum(value) - 0.5 * value)
  
  # Class
  cn = classification_normal(dataset$TIN)
  ct = classification_tumour(dataset$TIT)
  
  # Plot
  figure = ggplot(data = vtin,
         aes(x = '1', y = value, fill = variable)) +
    geom_bar(stat = 'identity') +
    mobster:::my_ggplot_theme() +
    facet_wrap(~ sample, ncol = 2, nrow = 1) +
    # coord_flip() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    coord_polar("y", start = 0, clip = 'off') +
    labs(
      title = bquote(bold('Sample composition')),
      x = element_blank(),
      y = "Percentage",
      subtitle = paste0(
        "TIN: ",
        round(vtin$value[1] * 100, 2),
        '%',
        " (",
        round(dataset$TIN_rf * 100, 2),
        '% RF) - ', cn['QC'],
        "\nTIT : ",
        round(vtin$value[3] * 100, 2),
        '% (',
        round(dataset$TIT_rf * 100, 2),
        '% RF) - ', ct['QC']
      )
    ) +
    scale_fill_manual(values = c(`Tumour cells` = 'plum4', `Normal cells` = 'plum2')) +
    guides(fill = guide_legend('')) +
    geom_label(
      aes(y = lab.ypos,
          label = paste0(round(value * 100, 1),  '%')),
      color = "black",
      fill = 'white'
      # fill = c(cn['color'], 'white', ct['color'], 'white')
    ) +
    theme(# strip.background = element_rect(
      #   color = 'black',
      #   fill = c(cn['color'], ct['color']),
      #   size=1.5,
      #   linetype="solid"
      # ),
      panel.background = element_rect(
        fill = c(cn['color'], ct['color']) %>% alpha(alpha = .4),
        # colour = cn['color'] %>% alpha(alpha = .4),
        size = 0.5,
        linetype = "solid"
      ))
  
  # 
  # normal = ggplot(vtin %>% filter(sample == 'Normal sample'),
  #                 aes(x = "", y = value, fill = variable)) +
  #   geom_bar(width = 1,
  #            stat = "identity",
  #            color = "white") +
  #   coord_polar("y", start = 0) +
  #   geom_label(aes(y = lab.ypos,
  #                  label = paste0(round(value * 100, 1),  '%')),
  #              color = "black",
  #              fill = 'white') +
  #   scale_fill_manual(values = c(`Tumour cells` = 'plum4', `Normal cells` = 'plum2')) +
  #   my_ggplot_theme() +
  #   guides(fill = guide_legend('')) +
  #   labs(
  #     x = '',
  #     y = '',
  #     title = bquote(bold(Normal)),
  #     subtitle = bquote(.(cn['QC']))
  #   ) +
  #   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  #   theme(
  #     legend.position = 'right',
  #     panel.background = element_rect(
  #       fill = cn['color'] %>% alpha(alpha = .4),
  #       colour = cn['color'] %>% alpha(alpha = .4),
  #       size = 0.5,
  #       linetype = "solid"
  #     ),
  #     panel.grid = element_blank(),
  #     plot.subtitle =  element_text(colour = cn['color'])
  #   )
  # 
  # tumour = ggplot(vtin %>% filter(sample == 'Tumour sample'),
  #                 aes(x = "", y = value, fill = variable)) +
  #   geom_bar(width = 1,
  #            stat = "identity",
  #            color = "white") +
  #   coord_polar("y", start = 0) +
  #   geom_label(aes(y = lab.ypos, label = paste0(round(value * 100, 1),  '%')),
  #              color = "black",
  #              fill = 'white') +
  #   scale_fill_manual(values = c(`Tumour cells` = 'plum4', `Normal cells` = 'plum2')) +
  #   my_ggplot_theme() +
  #   guides(fill = guide_legend('')) +
  #   labs(
  #     x = '',
  #     y = '',
  #     title = bquote(bold(Tumour)),
  #     subtitle = bquote(.(ct['QC']))
  #   ) +
  #   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  #   theme(
  #     legend.position = 'right',
  #     panel.background = element_rect(
  #       fill = ct['color'] %>% alpha(alpha = .4),
  #       colour = ct['color'] %>% alpha(alpha = .4),
  #       size = 0.5,
  #       linetype = "solid"
  #     ),
  #     panel.grid = element_blank(),
  #     plot.subtitle = element_text(colour = ct['color'])
  #   )
  # 
  # if (assemble)
  #   figure = ggpubr::ggarrange(
  #     tumour,
  #     normal,
  #     ncol = 1,
  #     nrow = 2,
  #     common.legend = TRUE,
  #     legend = 'bottom'
  #   )
  # else
  #   figure = list(normal = normal, tumour = tumour)
  
  return(figure)
  #
  #   cowplot::plot_grid(
  #     tumour,
  #     normal,
  #     ncol = 1,
  #     nrow = 2,
  #     align = 'v'
  #   )
  
}

classification_normal = function(TIN, console = FALSE)
{
  if (TIN < 0.01)
    return(c(
      `level` = 1,
      `color` = "forestgreen",
      `QC` = "No Contamination (<1%)"
    ))
  if (TIN < 0.03)
    return(c(
      `level` = 2,
      `color` = "steelblue",
      `QC` = "Minimal contamination (1-3%)"
    ))
  if (TIN < 0.07)
    return(c(
      `level` = 3,
      `color` = "goldenrod4",
      `QC` = "Some contamination (3-7%)"
    ))
  if (TIN < 0.15)
    return(c(
      `level` = 4,
      `color` = "indianred3",
      `QC` = "Contamination (7-15%)"
    ))
  return(c(
    `level` = 5,
    `color` = "purple",
    `QC` = "Huge contamination (>15%)"
  ))
}

classification_tumour = function(TIT)
{
  if (TIT < 0.15)
    return(c(
      `level` = 1,
      `color` = "purple",
      `QC` = "Impure (<15%)"
    ))
  if (TIT < 0.4)
    return(c(
      `level` = 1,
      `color` = "indianred3",
      `QC` = "Bad purity (15-45%)"
    ))
  if (TIT < 0.65)
    return(c(
      `level` = 3,
      `color` = "goldenrod1",
      `QC` = "Average purity (45-65%)"
    ))
  if (TIT < 0.85)
    return(c(
      `level` = 2,
      `color` = "steelblue",
      `QC` = "Good purity (65-85%)"
    ))
  return(c(
    `level` = 1,
    `color` = "forestgreen",
    `QC` = "High purity (>85%)"
  ))
}


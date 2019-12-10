plot_contamination_full_size = function(dataset)
{
  # fot the plots we half the estimates
  TIT = dataset$mobster_analysis$estimated_purity / 2
  TIN = dataset$BMix_analysis$estimated_purity / 2

  ggplot(
    dataset$joint %>% mutate(
      is_clonal = ifelse(is_clonal, 'Clonal', 'Other'),
      is_clonal = ifelse(VAF.tumour > 0, is_clonal, 'Excluded')
    ),
    aes(x = VAF.normal, y = VAF.tumour, color = is_clonal)
  ) +
    geom_point(
      data = dataset$joint,
      size = .3,
      alpha = 0.1,
      color = 'black'
    ) +
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
    mobster:::my_ggplot_theme() +
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
  # fot the plots we half the estimates
  TIT = dataset$mobster_analysis$estimated_purity / 2
  TIN = dataset$BMix_analysis$estimated_purity / 2


  data_table = dataset$joint %>%
    mutate(group = paste(tumour.cluster, '~', normal.cluster))

  group_labels = data_table$group %>% unique
  group_labels_colors = pio:::nmfy(group_labels, rep("black", length(group_labels)))

  for (cl in dataset$mobster_analysis$clonal_cluster)
    group_labels_colors[grepl(group_labels, pattern = cl)] = 'cyan4'


  ggplot(data_table,
         aes(
           x = VAF.normal,
           y = VAF.tumour,
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
    mobster:::my_ggplot_theme() +
    labs(title = paste0(""),
         x = 'VAF germline (zoom)',
         y = "VAF tumour (zoom)") +
    scale_color_manual(values = group_labels_colors)
}

plot_sample_contamination = function(dataset)
{
  # fot the plots we half the estimates
  TIT = dataset$mobster_analysis$estimated_purity / 2
  TIN = dataset$BMix_analysis$estimated_purity / 2

  vtin = data.frame(
    variable = c('Tumour cells', 'Normal cells', 'Tumour cells', 'Normal cells'),
    value = c(TIN * 2, 1 - (TIN * 2), TIT * 2, 1 - (TIT * 2)),
    sample = c(
      'Normal sample',
      'Normal sample',
      'Tumour sample',
      'Tumour sample'
    )
  )

  # ggplot(data = vtin,
  #        aes(x = '1', y = value, fill = variable)) +
  #   geom_bar(stat = 'identity') +
  #   mobster:::my_ggplot_theme() +
  #   facet_wrap(~ sample, ncol = 2, nrow = 1) +
  #   # coord_flip() +
  #   coord_polar("y", start = 0, clip = 'off') +
  #   labs(
  #     title = paste0("Sample composition"),
  #     x = element_blank(),
  #     y = "Percentage",
  #     subtitle = paste0(
  #       "TIN: ",
  #       round(vtin$value[1] * 100, 2),
  #       '%',
  #       " (",
  #       round(vtin$value[2] * 100, 2),
  #       '% germline)\n',
  #       "TIT : ",
  #       round(vtin$value[3] * 100, 2),
  #       '% (',
  #       round(vtin$value[4] * 100, 2),
  #       '% germline)'
  #     )
  #   ) +
  #   scale_fill_manual(values = c(`Tumour cells` = 'plum4', `Normal cells` = 'plum2')) +
  #   guides(fill = guide_legend(''))

  vtin = vtin %>%
    group_by(sample) %>%
    mutate(lab.ypos = cumsum(value) - 0.5 * value)

  cn = classification_normal(dataset$BMix_analysis$estimated_purity)
  ct = classification_tumour(dataset$mobster_analysis$estimated_purity)

  normal = ggplot(
    vtin %>% filter(sample == 'Normal sample'),
    aes(x = "", y = value, fill = variable)
  ) +
    geom_bar(width = 1,
             stat = "identity",
             color = "white") +
    coord_polar("y", start = 0) +
    geom_label(aes(y = lab.ypos, label = round(value * 100, 1)), color = "black", fill = 'white') +
    scale_fill_manual(values = c(`Tumour cells` = 'plum4', `Normal cells` = 'plum2')) +
    mobster:::my_ggplot_theme() +
    guides(fill = guide_legend('')) +
    labs(x = '', title = bquote(bold(Normal) ~ .(cn['QC']))) +
    theme(axis.text.x = NULL, axis.ticks.x = NULL) +
    theme(
      legend.position = 'right',
      panel.background = element_rect(
        fill = cn['color'] %>% alpha(alpha = .4),
        colour = cn['color'] %>% alpha(alpha = .4),
        size = 0.5,
        linetype = "solid"),
      plot.title = element_text(colour = cn['color'])
    )

  tumour = ggplot(
    vtin %>% filter(sample == 'Tumour sample'),
    aes(x = "", y = value, fill = variable)
  ) +
    geom_bar(width = 1,
             stat = "identity",
             color = "white") +
    coord_polar("y", start = 0) +
    geom_label(aes(y = lab.ypos, label = round(value * 100, 1)), color = "black", fill = 'white') +
    scale_fill_manual(values = c(`Tumour cells` = 'plum4', `Normal cells` = 'plum2')) +
    mobster:::my_ggplot_theme() +
    guides(fill = guide_legend('')) +
    labs(x = '', title = bquote(bold(Tumour) ~ .(ct['QC']))) +
    theme(axis.text.x = NULL, axis.ticks.x = NULL) +
    theme(
      legend.position = 'right',
      panel.background = element_rect(
        fill = ct['color'] %>% alpha(alpha = .4),
        colour = ct['color'] %>% alpha(alpha = .4),
        size = 0.5,
        linetype = "solid"),
      plot.title = element_text(colour = ct['color'])
    )


  cowplot::plot_grid(
    tumour,
    normal,
    ncol = 1,
    nrow = 2,
    align = 'v'
  )

}

classification_normal = function(TIN, console = FALSE)
{
  if(TIN < 0.01) return(c(`level` = 1, `color` = "forestgreen", `QC` = "No Contamination <1%"))
  if(TIN < 0.03) return(c(`level` = 2, `color` = "steelblue", `QC` = "Minimal contamination (1-3%)"))
  if(TIN < 0.07) return(c(`level` = 3, `color` = "goldenrod1", `QC` = "Some contamination (3-7%)"))
  if(TIN < 0.15) return(c(`level` = 4, `color` = "indianred3", `QC` = "Contamination (7-15%)"))
  return(c(`level` = 5, `color` = "purple", `QC` = "Huge contamination (>15%)"))
}

classification_tumour = function(TIT)
{
  if(TIT < 0.15) return(c(`level` = 1, `color` = "purple", `QC` = "Impure (<15%)"))
  if(TIT < 0.4) return(c(`level` = 1, `color` = "indianred3", `QC` = "Bad purity (15-45%)"))
  if(TIT < 0.65) return(c(`level` = 3, `color` = "goldenrod1", `QC` = "Average purity (45-65%)"))
  if(TIT < 0.85) return(c(`level` = 2, `color` = "steelblue", `QC` = "Good purity (65-85%)"))
  return(c(`level` = 1, `color` = "forestgreen", `QC` = "High purity (>85%)"))
}

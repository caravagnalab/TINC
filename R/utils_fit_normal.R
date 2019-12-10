plot_Bmix_fit = function(fit_normal, df_input, clonal_score){
  #  Plots
  nzero = sum(df_input$NV == 0)


  # plot_input = ggplot(x, aes(VAF)) + geom_histogram(binwidth = 0.01) + xlim(0, 1)
  plot_fit_normal = BMix::plot_clusters(fit_normal, data = df_input) +
    xlim(-0.01, 1) +
    labs(
      title = paste0("BMix fit - n = ", nrow(df_input), " (", nzero, ' zeros)'),
      subtitle = paste0(names(fit_normal$pi), ' (', round(fit_normal$pi, 2), '%)', collapse = ', ')
    ) +
    geom_vline(xintercept = BMix::Parameters(fit_normal) %>% pull(mean), linetype = 'dashed', size = .3) +
    ggrepel::geom_label_repel(
      data = BMix::Parameters(fit_normal),
      aes(x = mean, y = Inf, label = round(mean, 2)),
      fill = 'black', color = 'white', size = 3
    ) +
    scale_fill_brewer(palette = "Set2")

  # xlim(-0.01, min(1, max(df_input$NV/df_input$DP) * 2.5))

  # plot_fit_normal_d = BMix::plot_density(fit_normal, data = df_input)
  plot_fit_normal_d = bin_heatmap(C = median(df_input$DP), p = round(clonal_score, 2))

   cowplot::plot_grid(
    plot_fit_normal,
    plot_fit_normal_d,
    ncol = 2,
    nrow = 1,
    align = 'h'
  )
}


bin_heatmap = function(C = 42, offset = 1.5, p = 0.07)
{
  CMax = C * offset

  points = expand.grid(NV = 0:CMax, C = 1:CMax) %>%
    group_by(C) %>%
    filter(NV <= C/2) %>%
    data.frame

  points$z = apply(points, 1, function(x){
    dbinom(x[1], x[2], p)
  })

  max_C_plot = qbinom(p = .95, size = CMax, prob = p)

  q = ceiling(C*p)

  ggplot(points, aes(x = NV, y = C, fill = z)) +
    geom_raster() +
    coord_flip() +
    mobster:::my_ggplot_theme() +
    scale_x_continuous(limits = c(0, max_C_plot)) +
    # scale_fill_gradientn(
    #   colours = c("#0D0887FF", "#CC4678FF", "#F0F921FF")) +
    scale_fill_viridis_c(direction = -1, option = 'magma') +
    labs(
      y = "Depth at locus",
      x = "Number of reads with variant",
      title = 'Binomial density',
      subtitle = paste0("Coverage ", C, 'x (med.) ~ p = ', p)
    ) +
    guides(fill = guide_colorbar('Density', barwidth = 6)) +
    geom_hline(yintercept = C, size = .3, linetype = 'dashed') +
    geom_point(data = data.frame(x = q),
               inherit.aes = FALSE,
               aes(x = x, y = C),
               size = 4)




}

# TIN estimation
estimate_TIN = function(fit_normal, clonal_score)
{
  clonal_cluster = names(fit_normal$pi)
  # ccl =   paste0(round(fit_normal$pi, 2), ' x ', round(fit_normal$B.params, 2), collapse = ' + ')
  # clonal_cluster = names(which.max(fit_normal$B.params))
  # highest_Binomial_peak = fit_normal$B.params[clonal_cluster]
  highest_Binomial_peak = as.vector(clonal_score)

  # estimated_normal_purity = 1 - (0.5 - highest_Binomial_peak) * 2
  estimated_normal_purity =  highest_Binomial_peak * 2

  return(
    list(
      clonal_cluster = clonal_cluster,
      estimated_purity = estimated_normal_purity
    )
  )

}

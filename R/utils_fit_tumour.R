guess_mobster_clonal_cluster = function(mobster_fit_tumour, cutoff_miscalled_clonal)
{
  # Clonal cluster is guessed to be the highest peak below cutoff_miscalled_clonal = 510% VAF - seems reasonable without CNA
  all_clusters = mobster_fit_tumour$best$Clusters %>%
    filter(type == 'Mean', fit.value < cutoff_miscalled_clonal) %>%
    arrange(desc(fit.value)) %>%
    pull(cluster)

  # We also check that the clonal cluster has higher dimension of the others
  cl_map = sapply(all_clusters, mobster:::is_reasonable_clonal_cluster,
                  x = mobster_fit_tumour$best)

  clonal_cluster_fit = NULL
  if(all(!cl_map)) {
    warning("Check clonal cluster")
    clonal_cluster_fit = all_clusters[1]
  }
  else{
    if(!cl_map[1]) message("Not using C1 as clonal cluster")
    clonal_cluster_fit = names(which.max(cl_map))
  }

  return(clonal_cluster_fit)
  # # top-2 if Betas
  # if (length(all_clusters) > 1)
  # {
  #   top2 = all_clusters[1:2]
  #
  #   if (top2[2] != 'Tail')
  #   {
  #     m_top2 = mobster_fit_tumour$best$Clusters %>%
  #       filter(type == 'Mean', fit.value < cutoff_miscalled_clonal) %>%
  #       arrange(desc(fit.value))  %>%
  #       pull(fit.value)
  #
  #     d_m1 = mobster::ddbpmm(mobster_fit_tumour$best, data = m_top2[1])
  #     d_m2 = mobster::ddbpmm(mobster_fit_tumour$best, data = m_top2[2])
  #
  #     if (d_m2 < d_m1)
  #       clonal_cluster_fit = top2
  #   }
  # }
  #
  # clonal_cluster_fit = estimated_tumour_purity = clonal_tumour = NULL
  # m = 1
  #
  # repeat
  # {
  #   # Candidate clonal cluster
  #   clonal_cluster_fit = all_clusters[1:m]


  # # Stop with > 5 clonal mutations
  # if(length(clonal_tumour) > 5) break;
  #
  # Stop if we tried all
  #   if(m == length(all_clusters)) stop("No clonal mutations with this LV value, try lowering.")
  #
  #   message("Augmenting clonal cluster to achieve more mutations")
  #   m = m + 1
  # }
}

plot_mobster_fit = function(mobster_fit_tumour, cutoff_lv_assignment, clonal_cluster)
{
  # plot_input = ggplot(x, aes(VAF)) + geom_histogram(binwidth = 0.01) + xlim(0, 1)
  plot_fit_tumour = mobster::plot.dbpmm(mobster_fit_tumour$best, cutoff_assignment = cutoff_lv_assignment) +
    labs(title = paste0('MOBSTER (clonal: ', clonal_cluster, ')'), caption = NULL)

  plot_lv = mobster::plot_latent_variables(mobster_fit_tumour$best, cutoff_lv_assignment)
  # plot_fs = mobster::plot_entropy(mobster_fit_tumour$best)

  cowplot::plot_grid(
    # plot_input,
    plot_fit_tumour,
    plot_lv,
    # plot_fs,
    ncol = 2,
    nrow = 1,
    align = 'h',
    axis = 'bt'
  )
}

# The clonal cluster is
# - the highest-mean peak below cutoff_miscalled_clonal;
# - due to a Beta distribution;
# - non-gathered on chromosome locationns, if geonme locations are available
guess_mobster_clonal_cluster = function(mobster_fit_tumour,
                                        cutoff_miscalled_clonal,
                                        chromosomes = paste0('chr', 1:22))
{
  # Clonal cluster is guessed to be the highest peak below cutoff_miscalled_clonal, which- seems reasonable without CNA
  all_clusters = mobster_fit_tumour$best$Clusters %>%
    filter(type == 'Mean', cluster != 'Tail', fit.value < cutoff_miscalled_clonal) %>%
    arrange(desc(fit.value))

  # Want them to be > 0.05 to denote the fact that it is uniform across `chromosomes`
  all_clusters$chisq_pval = sapply(
    all_clusters %>%
      pull(cluster),
    location_likelihood,
    mobster_fit_tumour = mobster_fit_tumour,
    chromosomes = chromosomes)

  all_clusters$OK_chisq = all_clusters$chisq_pval > 0.05

  # If there is one cluster -- return it whatever it is
  if(length(all_clusters) == 1) {
    if(!all_clusters$OK_chisq[1]) {
      cli::cli_alert_danger("One clonal cluster, failing the Chi square test for genomic locations")
    }
    else
      cli::cli_alert_success("One clonal cluster, passing the Chi square test for genomic locations")

    return(all_clusters$cluster)
  }

  # If there are more, and all non OK_chisq, return top
  if(all(!all_clusters$OK_chisq)) {
    cli::cli_alert_danger("All clusters of mutations fail the Chi square test for genomic locations")

    return(all_clusters$cluster[1])
  }

  # Returnn top OK_chisq
  rank = all_clusters %>%
    dplyr::filter(OK_chisq) %>%
    dplyr::pull(cluster)

  return(rank)

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

location_likelihood = function(mobster_fit_tumour, cluster, chromosomes = paste0('chr', 1:22))
{
  assignments = mobster_fit_tumour$best$data %>%
    dplyr::filter(cluster == !!cluster)

  assignments = CNAqc:::relative_to_absolute_coordinates(assignments)

  coordinates = CNAqc::chr_coordinates_hg19 %>%
    dplyr::filter(chr %in% chromosomes)

  # N outcomes ~ one per chromosome
  n = nrow(coordinates)

  all_mapping = rep(0, n)
  names(all_mapping) = paste0(coordinates$chr)

  # Observations (counts)
  mapping_counts = table(assignments$chr)
  all_mapping[names(mapping_counts)] = mapping_counts

  # Chisquare assuming p is uniform
  test = chisq.test(all_mapping, p = rep(1/n, n))

  test$p.value
  # n = nrow(CNAqc::chr_coordinates_hg19) * 2

  # names(all_mapping) = c(
  #   paste0(CNAqc::chr_coordinates_hg19$chr, 'p'),
  #   paste0(CNAqc::chr_coordinates_hg19$chr, 'q')
  # )

  # st_chr = pio:::nmfy(
  #   CNAqc::chr_coordinates_hg19 %>%pull(chr),
  #   CNAqc::chr_coordinates_hg19 %>%pull(centromerStart)
  # )
  #
  # stp_chr = pio:::nmfy(
  #   CNAqc::chr_coordinates_hg19 %>%pull(chr),
  #   CNAqc::chr_coordinates_hg19 %>%pull(centromerEnd)
  # )

  # mapping = assignments %>%
  #   mutate(
  #     arm = ifelse(
  #       from < st_chr[chr] & from < stp_chr[chr],
  #       "p",
  #       "q"
  #     ),
  #     label = paste0(chr, arm)
  #   ) %>%
  #   pull(label)

  # mapping_counts = table(mapping)
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

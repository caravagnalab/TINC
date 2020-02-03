# The clonal cluster is
# - the highest-mean peak below cutoff_miscalled_clonal;
# - due to a Beta distribution;
# - non-gathered on chromosome locationns, if geonme locations are available
guess_mobster_clonal_cluster = function(mobster_fit_tumour,
                                        cutoff_miscalled_clonal,
                                        chromosomes = paste0('chr', 1:22))
{
  # Clonal cluster is guessed to be the highest non-tail peak below cutoff_miscalled_clonal, which - seems reasonable without CNA
  all_clusters = mobster_fit_tumour$best$Clusters %>%
    filter(type == 'Mean', cluster != 'Tail', fit.value < cutoff_miscalled_clonal) %>%
    arrange(desc(fit.value))

  # Handle the extreme case of only one Beta cluster above cutoff_miscalled_clonal
  if(length(all_clusters) == 0) {
    # Force picking that cluster anyway
    all_clusters = mobster_fit_tumour$best$Clusters %>%
      filter(type == 'Mean', cluster != 'Tail') %>%
      arrange(desc(fit.value))
  }

  cat('\n')

  # Apply the location_likelihood_test to check that mutationns spread evenly §across `chromosomes`
  all_clusters$location_likelihood_test = sapply(
    all_clusters %>%
      pull(cluster),
    location_likelihood,
    mobster_fit_tumour = mobster_fit_tumour,
    chromosomes = chromosomes)

  # all_clusters$OK_chisq = all_clusters$chisq_pval > 0.05
  all_clusters$OK_8020 = all_clusters$location_likelihood_test

  # If there is one cluster -- return it whatever it is
  if(length(all_clusters) == 1) {
    if(!all_clusters$OK_8020[1]) {
      cli::cli_alert_danger("Location Likelihood: one clonal cluster that fails the test, non better choice.")
    }
    else
      cli::cli_alert_success("Location Likelihood: one clonal cluster, passing the test.")

    return(all_clusters$cluster)
  }

  # If there are more, and all non OK_chisq, return top
  if(all(!all_clusters$OK_8020)) {
    cli::cli_alert_danger("Location Likelihood: all clusters  fail the test, returning the first one.")

    return(all_clusters$cluster[1])
  }

  # Return top OK
  rank = all_clusters %>%
    dplyr::filter(OK_8020) %>%
    dplyr::pull(cluster)

  if(rank != "C1")
    cli::cli_alert_danger("Location Likelihood: changed C1 to cluster {.value {rank[1]}}")

  return(rank[1])

  # We also check that the clonal cluster has higher dimension of the others
  # cl_map = sapply(all_clusters, mobster:::is_reasonable_clonal_cluster,
  #                 x = mobster_fit_tumour$best)
  #
  # clonal_cluster_fit = NULL
  # if(all(!cl_map)) {
  #   warning("Check clonal cluster")
  #   clonal_cluster_fit = all_clusters[1]
  # }
  # else{
  #   if(!cl_map[1]) message("Not using C1 as clonal cluster")
  #   clonal_cluster_fit = names(which.max(cl_map))
  # }
  #
  # return(clonal_cluster_fit)
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

# location_likelihood = function(mobster_fit_tumour, cluster, chromosomes)
# {
#   # Clusters
#   assignments = mobster_fit_tumour$best$data %>%
#     dplyr::filter(cluster == !!cluster)
#
#   # Chromosomes used (actual mutations + a list)
#   used_chromosomes = unique(c(assignments$chr, chromosomes))
#
#   assignments = CNAqc:::relative_to_absolute_coordinates(assignments)
#
#   # Coordinates used
#   coordinates = CNAqc::chr_coordinates_hg19 %>%
#       dplyr::filter(chr %in% used_chromosomes)
#
#   G = sum(coordinates$length)
#   coordinates$p = coordinates$length/G
#
#   # N outcomes ~ one per chromosome
#   n = nrow(coordinates)
#
#   all_mapping = rep(0, n)
#   names(all_mapping) = paste0(coordinates$chr)
#
#   # Observations (counts)
#   mapping_counts = table(assignments$chr)
#   all_mapping[names(mapping_counts)] = mapping_counts
#
#   # Chisquare assuming p is uniform
#   test = chisq.test(all_mapping, p = coordinates$p)
#
#   cli::cli_alert_info(
#     "Cluster {.value {cluster}} has Chi-square p-value p = {.value {test$p.value}}"
#   )
#   print(all_mapping)
#
#   test$p.value
#   # n = nrow(CNAqc::chr_coordinates_hg19) * 2
#
#   # names(all_mapping) = c(
#   #   paste0(CNAqc::chr_coordinates_hg19$chr, 'p'),
#   #   paste0(CNAqc::chr_coordinates_hg19$chr, 'q')
#   # )
#
#   # st_chr = pio:::nmfy(
#   #   CNAqc::chr_coordinates_hg19 %>%pull(chr),
#   #   CNAqc::chr_coordinates_hg19 %>%pull(centromerStart)
#   # )
#   #
#   # stp_chr = pio:::nmfy(
#   #   CNAqc::chr_coordinates_hg19 %>%pull(chr),
#   #   CNAqc::chr_coordinates_hg19 %>%pull(centromerEnd)
#   # )
#
#   # mapping = assignments %>%
#   #   mutate(
#   #     arm = ifelse(
#   #       from < st_chr[chr] & from < stp_chr[chr],
#   #       "p",
#   #       "q"
#   #     ),
#   #     label = paste0(chr, arm)
#   #   ) %>%
#   #   pull(label)
#
#   # mapping_counts = table(mapping)
# }

location_likelihood = function(mobster_fit_tumour, cluster, chromosomes)
{
  # Clusters
  assignments = mobster_fit_tumour$best$data %>%
    dplyr::filter(cluster == !!cluster)

  # Chromosomes used (actual mutations + a list)
  used_chromosomes = unique(c(assignments$chr, chromosomes))

  # Counts for all used chromosomes
  all_mapping = rep(0, length(used_chromosomes))
  names(all_mapping) = used_chromosomes

  mapping_counts = table(assignments$chr)
  all_mapping[names(mapping_counts)] = mapping_counts

  # Sort them by size, and determinne 80% prob
  all_mapping = sort(all_mapping, decreasing = TRUE)

  min_cutoff_mapping = sum(all_mapping) * 0.01
  if(min_cutoff_mapping > 10) min_cutoff_mapping = 10

  if(any(all_mapping <= min_cutoff_mapping, na.rm = T))
  {
    cli::cli_alert_warning("Location likelihood: chromosomes arms have < {.value {min_cutoff_mapping + 1}} mutations, and will not be considered.")
    print(all_mapping)

    all_mapping = all_mapping[all_mapping > min_cutoff_mapping]
  }

  total_mass = sum(all_mapping)
  cumulative_mass = cumsum(all_mapping)


  index_80p = sum(cumulative_mass < total_mass * 0.6) + 1
  index_20p = length(used_chromosomes) * .2

  if(index_80p > index_20p)
    cli::cli_alert_success(
      "Cluster {.value {cluster}}: 60% of counts in {.value {index_80p}} chromosomes (spread over >20% of used chromosomes: {.value {index_20p}})."
    )
  else
    cli::cli_alert_warning(
      "Cluster {.value {cluster}}: 60% of counts in {.value {index_80p}} chromosomes (spread in less than20% of used chromosomes: {.value {index_20p}}) -- is this cluster do to CNA?"
    )

  # TRUE if it is OK
  return(index_80p > index_20p)
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

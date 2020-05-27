analysis_mode = function(x)
{
  if(all(is.null(x))) return("NO_CNA")

  return("CNA")
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Analyse tumour sample with MOBSTER
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
analyse_mobster = function(x,
                           cutoff_miscalled_clonal,
                           cna_map,
                           cutoff_lv_assignment,
                           chromosomes,
                           ...)
{
  cli::cli_h1("Analysing tumour sample with MOBSTER")
  cat('\n')

  #  MOBSTER fit of the input data
  mobster_fit_tumour = mobster::mobster_fit(x,
                                            ...)

  output = mobster::Clusters(mobster_fit_tumour$best) %>%
    rename(tumour.cluster = cluster)

  # Determine clonal cluster
  clonal_cluster = TINC:::guess_mobster_clonal_cluster(
    mobster_fit_tumour = mobster_fit_tumour,
    cutoff_miscalled_clonal = cutoff_miscalled_clonal,
    use_heuristic = TINC:::analysis_mode(cna_map) == "NO_CNA",
    chromosomes = chromosomes)

  stopifnot(length(clonal_cluster) > 0)

  ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
  # Tumour purity: estimated based on the fact that we could map or no the mutations to CNA
  estimated_tumour_purity = params_purity_adj = NULL

  if(all(is.null(cna_map)))
  {
    cli::cli_alert_info("Without CNA, TINC will be estimating tumour purity as 2*x, with x the clonal peak.")

    # If we could not, then we just assume everything is diploid, therefore
    # therefore 2 * the Beta peak of the clonal cluster

    ccluster_mean = mobster_fit_tumour$best$Clusters %>%
      filter(cluster %in% clonal_cluster, type == 'Mean') %>%
      pull(fit.value) %>% mean

    estimated_tumour_purity = ccluster_mean * 2

    params_purity_adj = list(
      mean_cluster_VAF = ccluster_mean, # Observed VAF
      karyotype = NA,
      minor_copies = NA,
      Major_copies = NA,
      mutant_alleles = NA
    )
  }
  else
  {
    # Otherwise we do something a little bit smarter, which is normalising for segment's ploidy.
    cli::cli_alert_info("With CNA, TINC will be wstimating tumour purity as 2*x, with x the clonal peak.")

    # We take the mean of the clonal cluster
    ccluster_mean = mobster_fit_tumour$best$Clusters %>%
      filter(cluster %in% clonal_cluster, type == 'Mean') %>%
      pull(fit.value) %>%
      mean

    # We found the karyotype of the mutations that we used (we take the first one, because they are all the same)
    used_karyotype = strsplit(x$karyotype[1], ':')[[1]]

    # By design, we should have correctly identified the clonal cluster, as the set of
    # mutations that happened before aneuploidy. Also, again by the fact that we use
    # only simple karyotypes, we assume that the mutation is present in a number of
    # copies that match the actual Major allele (M).
    mut.allele = case_when(
      x$karyotype[1] == "1:0" ~ 1,
      x$karyotype[1] == "2:0" ~ 2,
      x$karyotype[1] == "1:1" ~ 1,
      x$karyotype[1] == "2:1" ~ 2,
      x$karyotype[1] == "2:2" ~ 2,
      TRUE ~ -1
     )

    if(mut.allele == -1) stop("Karyotype not recognised: ", x$karyotype[1])

    estimated_tumour_purity = CNAqc:::purity_estimation_fun(
      v = ccluster_mean, # Observed VAF
      m = used_karyotype[2] %>% as.numeric,
      M = used_karyotype[1] %>% as.numeric,
      mut.allele = mut.allele
      )

    params_purity_adj = list(
      mean_cluster_VAF = ccluster_mean, # Observed VAF
      karyotype = paste(used_karyotype, collapse =':'),
      minor_copies = used_karyotype[2] %>% as.numeric,
      Major_copies = used_karyotype[1] %>% as.numeric,
      mutant_alleles = mut.allele
    )
  }
  ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

  # List of clonal mutations in the tumour, with LV > cutoff_lv_assignment
  clonal_tumour = NULL
  repeat
  {
    clonal_tumour = mobster::Clusters(mobster_fit_tumour$best, cutoff_assignment = cutoff_lv_assignment) %>%
      filter(cluster %in% clonal_cluster) %>%
      pull(id)

    if(length(clonal_tumour) > 20 | cutoff_lv_assignment < 0) break

     cutoff_lv_assignment = cutoff_lv_assignment - 0.03

     # cat("Dynamic adjustment: ", cutoff_lv_assignment, " n = ", length(clonal_tumour))
  }

  # Plot
  figure = plot_mobster_fit(mobster_fit_tumour, cutoff_lv_assignment, clonal_cluster)

  return(
    list(
      output = output,
      fit = mobster_fit_tumour,
      plot = figure,
      clonal_cluster = clonal_cluster,
      clonal_mutations = clonal_tumour,
      estimated_purity = estimated_tumour_purity,
      params_purity_adj = params_purity_adj,
      cutoff_lv_assignment = cutoff_lv_assignment
    )
  )
}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Analyse normal samples with BMix
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
analyse_BMix = function(
  x,
  cna_map,
  ...)
{
  cli::cli_h1("Analysing normal sample with BMix")
  cat('\n')

  if(nrow(x) == 0)
  {
    stop("There are no tumour clonal mutations in the normal sample, there is no contamination?")
  }

  df_input = x %>%
    dplyr::select(NV, DP) %>%
    as.data.frame

  fit_normal = BMix::bmixfit(df_input, ...)

  output = x
  output$normal.cluster = fit_normal$labels

  Binomial_peaks = BMix::Parameters(fit_normal) %>% pull(mean)
  Binomial_pi = fit_normal$pi

  clonal_score = (Binomial_peaks %*% Binomial_pi) %>% as.numeric

  # options(warn=-1)
  figure = plot_Bmix_fit(fit_normal, df_input, clonal_score)
  # options(warn=0)

  # TIN = estimate_TIN(fit_normal, clonal_score, cna_map)
  #
  clonal_cluster = names(fit_normal$pi)
  highest_Binomial_peak = as.vector(clonal_score)

  ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######
  # TIN Contamination is a form of purity: it is estimated as for cancer
  estimated_normal_purity = params_purity_adj = NULL

  if(all(is.null(cna_map)))
  {
    # Without CNA it is
    estimated_normal_purity =  highest_Binomial_peak * 2

    params_purity_adj = list(
      mean_cluster_VAF = highest_Binomial_peak, # Observed VAF
      karyotype = NA,
      minor_copies = NA,
      Major_copies = NA,
      mutant_alleles = NA
    )
  }
  else
  {
    # Otherwise we use CNA as for tumour purity
    used_karyotype = strsplit(x$karyotype[1], ':')[[1]]

    # We behave as for mobster fits.
    mut.allele = case_when(
      x$karyotype[1] == "1:0" ~ 1,
      x$karyotype[1] == "2:0" ~ 2,
      x$karyotype[1] == "1:1" ~ 1,
      x$karyotype[1] == "2:1" ~ 2,
      x$karyotype[1] == "2:2" ~ 2,
      TRUE ~ -1
    )

    if(mut.allele == -1) stop("Karyotype not recognised: ", x$karyotype[1])

    estimated_normal_purity = CNAqc:::purity_estimation_fun(
      v = highest_Binomial_peak, # Observed VAF
      m = used_karyotype[2] %>% as.numeric,
      M = used_karyotype[1] %>% as.numeric,
      mut.allele = mut.allele
    )

    params_purity_adj = list(
      mean_cluster_VAF = highest_Binomial_peak, # Observed VAF
      karyotype = paste(used_karyotype, collapse =':'),
      minor_copies = used_karyotype[2] %>% as.numeric,
      Major_copies = used_karyotype[1] %>% as.numeric,
      mutant_alleles = mut.allele
    )
  }
  ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### #######

  #
  clonal_mutations = output %>%
    filter(normal.cluster %in% clonal_cluster) %>% pull(id)

  # pio::pioStr("\n    Binomial peaks", paste(Binomial_peaks, collapse = ' '), suffix = '\n')
  # pio::pioStr("Mixing proportions", paste(Binomial_pi, collapse = ' '), suffix = '\n')
  # pio::pioStr("      Clonal score", clonal_score, suffix = '\n')
  # pio::pioStr("               TIN", TIN$estimated_purity, suffix = '\n')

  cli::cli_alert_success(
    "Binomial peaks {.value {Binomial_peaks}} with proportions {.value {Binomial_pi}}. Clonal score {.value {clonal_score}} with TINN {.value {estimated_normal_purity}}"
  )

  return(
    list(
      output = output,
      fit = fit_normal,
      plot = figure,
      clonal_cluster = clonal_cluster,
      clonal_mutations = clonal_mutations,
      estimated_purity = estimated_normal_purity,
      params_purity_adj = params_purity_adj
    )
  )

}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Analyse tumour sample with VIBER
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
analyze_VIBER = function(x, ...)
{
  cli::cli_h1("Analysing tumour and normal samples with VIBER")
  cat('\n')

  options(easypar.parallel = FALSE)

  # Only OK tumour
  NV.n = as_normal(x) %>% dplyr::filter(OK_tumour) %>% dplyr::select(NV) %>% rename(Normal = NV)
  NV.t = as_tumour(x) %>% dplyr::filter(OK_tumour) %>% dplyr::select(NV) %>% rename(Tumour = NV)

  DP.n = as_normal(x) %>% dplyr::filter(OK_tumour) %>% dplyr::select(DP) %>% rename(Normal = DP)
  DP.t = as_tumour(x) %>% dplyr::filter(OK_tumour) %>% dplyr::select(DP) %>% rename(Tumour = DP)

  fit = VIBER::variational_fit(
    dplyr::bind_cols(NV.n, NV.t),
    dplyr::bind_cols(DP.n, DP.t),
    ...
  )

  options(easypar.parallel = TRUE)

  cat('\n')

  fit = VIBER::choose_clusters(fit, binomial_cutoff = 0, dimensions_cutoff = 0, pi_cutoff = 0.02)
  pl = VIBER::plot_2D(fit, "Normal", "Tumour") + labs(title = 'VIBER analysis')

  # Table
  # x %>%
  #   full_join(
  #     dplyr::bind_cols(
  #       fit$labels %>% rename(VIBER.cluster = cluster.Binomial),
  #       x %>% dplyr::filter(OK_tumour) %>% dplyr::select(id)
  #     ),
  #     by = 'id')
  out = dplyr::bind_cols(
        fit$labels %>% rename(VIBER.cluster = cluster.Binomial),
        x %>% dplyr::filter(OK_tumour) %>% dplyr::select(id)
        )




  return(list(output = out, fit=fit, plot=pl))
}

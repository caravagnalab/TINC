# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Analyse tumour sample with MOBSTER
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
analyse_mobster = function(x,
                           cutoff_miscalled_clonal,
                           cutoff_lv_assignment,
                           chromosomes,
                           ...)
{
  cli::cli_h1("Analysing tumour sample with MOBSTER")
  cat('\n')

  #  MOBSTER fit
  mobster_fit_tumour = mobster::mobster_fit(x,
                                            ...)

  output = mobster::Clusters(mobster_fit_tumour$best) %>%
    rename(tumour.cluster = cluster)

  # Determine clonal cluster
  clonal_cluster = guess_mobster_clonal_cluster(mobster_fit_tumour = mobster_fit_tumour,
                                                cutoff_miscalled_clonal = cutoff_miscalled_clonal,
                                                chromosomes = chromosomes)

  stopifnot(length(clonal_cluster) > 0)


  # Tumour purity is therefore 2 * the Beta peak of the clonal cluster
  estimated_tumour_purity = mobster_fit_tumour$best$Clusters %>%
    filter(cluster %in% clonal_cluster, type == 'Mean') %>%
    pull(fit.value) %>% mean * 2

  # List of clonal mutations in the tumour, with LV > cutoff_lv_assignment
  clonal_tumour = mobster::Clusters(mobster_fit_tumour$best, cutoff_assignment = cutoff_lv_assignment) %>%
    filter(cluster %in% clonal_cluster) %>%
    pull(id)

  # Plot
  figure = plot_mobster_fit(mobster_fit_tumour, cutoff_lv_assignment, clonal_cluster)

  return(
    list(
      output = output,
      fit = mobster_fit_tumour,
      plot = figure,
      clonal_cluster = clonal_cluster,
      clonal_mutations = clonal_tumour,
      estimated_purity = estimated_tumour_purity
    )
  )
}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Analyse normal samples with BMix
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
analyse_BMix = function(x, ...)
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

  TIN = estimate_TIN(fit_normal, clonal_score)

  clonal_mutations = output %>%
    filter(normal.cluster %in% TIN$clonal_cluster) %>% pull(id)

  # pio::pioStr("\n    Binomial peaks", paste(Binomial_peaks, collapse = ' '), suffix = '\n')
  # pio::pioStr("Mixing proportions", paste(Binomial_pi, collapse = ' '), suffix = '\n')
  # pio::pioStr("      Clonal score", clonal_score, suffix = '\n')
  # pio::pioStr("               TIN", TIN$estimated_purity, suffix = '\n')

  cli::cli_alert_success(
    "Binomial peaks {.value {Binomial_peaks}} with proportions {.value {Binomial_pi}}. Clonal score {.value {clonal_score}} with TINN {.value {TIN$estimated_purity}}"
  )

  return(
    list(
      output = output,
      fit = fit_normal,
      plot = figure,
      clonal_cluster = TIN$clonal_cluster,
      clonal_mutations = clonal_mutations,
      estimated_purity = TIN$estimated_purity
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

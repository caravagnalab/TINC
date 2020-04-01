to_string = function(x)
{
  x$params$VAF_range_tumour = paste(x$params$VAF_range_tumour, collapse = ':')
  dfp = data.frame(x$params, stringsAsFactors = FALSE)
  colnames(dfp) = paste0('params_',   colnames(dfp))

  with_cna = !(all(is.null(x$fit$CNA)))

  dfp2 = NULL
  if(with_cna)
  {
    dfp2 = data.frame(x$cna_params, stringsAsFactors = FALSE)
    colnames(dfp2) = paste0('cna_params_',   colnames(dfp2))
  }

  data.frame(
    TIN = x$TIN,
    TIT = x$TIT,
    n_mutations_total = x$data %>% nrow,
    n_mutations_clonal = x$fit$mobster_analysis$clonal_mutations %>% length,
    clonal_cluster = x$fit$mobster_analysis$clonal_cluster,
    latentvar_clonal_cutoff = x$fit$mobster_analysis$cutoff_lv_assignment,
    n_binomial_clusters_bmix = x$fit$BMix_analysis$clonal_cluster %>% length,
    n_binomial_clusters_viber = x$fit$VIBER_analysis$fit$K,
    binomial_peaks_bmix = x$fit$BMix_analysis$fit$B.params %>% paste(collapse = ':'),
    stringsAsFactors = FALSE
  ) %>%
    bind_cols(dfp) %>%
    bind_cols(dfp2)
}

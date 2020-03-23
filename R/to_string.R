to_string = function(x)
{
  x$params$VAF_range_tumour = paste(x$params$VAF_range_tumour, collapse = ':')
  dfp = data.frame(x$params, stringsAsFactors = FALSE)
  colnames(dfp) = paste0('params_',   colnames(dfp))

  data.frame(
    TIN = x$TIN,
    TIT = x$TIT,
    n_mutations_total = x$data %>% nrow,
    n_mutations_clonal = x$fit$mobster_analysis$clonal_mutations %>% length,
    clonal_cluster = x$fit$mobster_analysis$clonal_cluster,
    latentvar_clonal_cutoff = x$fit$mobster_analysis$cutoff_lv_assignment,
    n_binomial_clusters_bmix = x$fit$BMix_analysis$clonal_cluster %>% length,
    n_binomial_clusters_viber = x$fit$VIBER_analysis$fit$K,
    stringsAsFactors = FALSE
  ) %>%
    bind_cols(dfp)
}

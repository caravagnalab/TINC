#' Create a data frame representation of a TINC fit
#'
#' @param x A TINC fit.
#'
#' @return A dataframe.
#' @export
#'
#' @examples
to_string = function(x)
{
  x$params$VAF_range_tumour = paste(x$params$VAF_range_tumour, collapse = ':')
  dfp = data.frame(x$params, stringsAsFactors = FALSE)
  colnames(dfp) = paste0('params_',   colnames(dfp))

  # With CNA data (augmented properly)
  with_cna = !(all(is.null(x$fit$CNA)))

  dfp2 = NULL
  if(with_cna)
  {
    dfp2 = data.frame(x$cna_params, stringsAsFactors = FALSE)
    colnames(dfp2) = paste0('cna_params_',   colnames(dfp2))
  }
  else
  {
    vv = rep(NA, 3)
    names(vv) = c("cna_params_used_chromosomes", "cna_params_karyotype", "cna_params_ploidy")
    dfp2 = t(data.frame(vv, stringsAsFactors = FALSE))
    dfp2 = data.frame(dfp2)
  }

  data.frame(
    TIN = x$TIN,
    TIT = x$TIT,
    TIN_rf = x$TIN_rf,
    TIT_rf = x$TIT_rf,
    n_mutations_total = x$data %>% nrow,
    n_mutations_clonal = x$fit$mobster_analysis$clonal_mutations %>% length,
    clonal_cluster = x$fit$mobster_analysis$clonal_cluster,
    latentvar_clonal_cutoff = x$fit$mobster_analysis$cutoff_lv_assignment,
    betas_mobster = x$fit$mobster_analysis$fit$best$Kbeta,
    tail_mobster = x$fit$mobster_analysis$fit$best$Kbeta,
    K_mobster = x$fit$mobster_analysis$fit$best$K,
    n_binomial_clusters_bmix = x$fit$BMix_analysis$clonal_cluster %>% length,
    n_binomial_clusters_viber = x$fit$VIBER_analysis$fit$K,
    binomial_peaks_bmix = x$fit$BMix_analysis$fit$B.params %>% paste(collapse = ':'),
    stringsAsFactors = FALSE
  ) %>%
    bind_cols(dfp) %>%
    bind_cols(dfp2)
}

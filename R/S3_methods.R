#' Print a TINC object to screen
#'
#' @param x A TINC analysis computed with \code{autofit}.
#' @param ... Extra S3 parameters
#'
#' @return
#' @export
#'
#' @examples
#' # Automatic call
#' autofit(random_TIN(), FAST = TRUE)
print.tin_obj = function(x, ...)
{
  if(!inherits(x, "tin_obj")) stop("Not a TINC object .... run autofit(.) first, aborting.")

  verbose = FALSE
  pio::pioHdr("TINC - Profiler for bulk samples TIN/TIT contamination")

  cat('\n')
  pio::pioStr(' Input : ',
              'n =', sum(x$fit$joint$used), 'used out of', nrow(x$fit$joint),
              paste0('annotated (', round(sum(x$fit$joint$used)/nrow(x$fit$joint), 2) * 100, '%)'), suffix = '\n')

  pio::pioStr('        TIT : ',
              paste0(round(x$fit$mobster_analysis$estimated_purity, 2) * 100, '%'),
              ' ~ n =', length(x$fit$mobster_analysis$clonal_mutations),
              'clonal mutations, cluster', x$fit$mobster_analysis$clonal_cluster,
              suffix = '\n')

  pio::pioStr('        TIN : ',
              paste0(round(x$fit$BMix_analysis$estimated_purity, 2) * 100, '%'),
              ' ~ n =',
              sum(x$fit$BMix_analysis$output$VAF > 0), 'with VAF > 0',
              suffix = '\n')

  # Report Classification
  cat('\n')

  cn = classification_normal(x$fit$BMix_analysis$estimated_purity)
  ct = classification_tumour(x$fit$mobster_analysis$estimated_purity)

  # pc = c('forestgreen', 'steelblue', 'goldenrod1', 'indianred3', 'purple')

  if(ct['color'] == 'forestgreen') cat(black$bgGreen$bold("   QC Tumour  "), ct['QC'] %>% green)
  if(ct['color'] == 'steelblue') cat(black$bgBlue$bold("   QC Tumour  "), ct['QC'] %>% blue)
  if(ct['color'] == 'goldenrod1') cat(black$bgYellow$bold("   QC Tumour  "), ct['QC'] %>% yellow)
  if(ct['color'] == 'indianred3') cat(black$bgRed$bold("   QC Tumour  "), ct['QC'] %>% red)
  if(ct['color'] == 'purple') cat(black$bgMagenta$bold("   QC Tumour  "), ct['QC'] %>% magenta)

  cat('\n')

  if(cn['color'] == 'forestgreen') cat(black$bgGreen$bold("   QC Normal  "), cn['QC'] %>% green)
  if(cn['color'] == 'steelblue') cat(black$bgBlue$bold("   QC Normal  "), cn['QC'] %>% blue)
  if(cn['color'] == 'goldenrod1') cat(black$bgYellow$bold("   QC Normal  "), cn['QC'] %>% yellow)
  if(cn['color'] == 'indianred3') cat(black$bgRed$bold("   QC Normal  "), cn['QC'] %>% red)
  if(cn['color'] == 'purple') cat(black$bgMagenta$bold("   QC Normal  "), cn['QC'] %>% magenta)

  if(verbose)
  {
    pio::pioTit('Profiled data:', x$fit$file)
    pio::pioDisp(x$fit$joint)

    pio::pioTit('MOBSTER analysis of tumour sample')
    print(x$fit$mobster_analysis$fit$best)

    pio::pioTit('BMix analysis of tumour sample')
    print(x$fit$BMix_analysis$fit)

    pio::pioTit('VIBER analysis of tumour sample')
    print(x$fit$VIBER_analysis$fit)
  }
}

#' Plot a TINC analysis.
#'
#' @param x A TINC analysis computed with \code{autofit}.
#' @param ... Extra S3 parameters.
#'
#' @return A TINC analysis plot (a `ggplot` figure with multiple panels).
#' @export
#'
#' @examples
#' plot(autofit(random_TIN(), FAST = TRUE))
plot.tin_obj = function(x, ...)
{
  if(!inherits(x, "tin_obj")) stop("Not a TINC object .... run autofit(.) first, aborting.")

  x$plot
}

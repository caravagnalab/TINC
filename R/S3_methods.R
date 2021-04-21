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
#' rt = random_TIN()
#' autofit(input = rt$data, cna = rt$cna, FAST = TRUE)
print.tin_obj = function(x, ...)
{
  if(!inherits(x, "tin_obj")) stop("Not a TINC object .... run autofit(.) first, aborting.")

  verbose = FALSE
  cli::cli_h1("TINC profiler for bulk WGS")
  cat('\n')

  if(!(all(is.null(x$fit$CNA))))
  {
    cli::cli_alert_info("Copy Number data has been used for this analysis (karyotype {.field {x$data$karyotype[1]}})")
    cat("\n")
    
    cli::cli_rule()
    x$fit$CNA %>% print
    cli::cli_rule()
  }
  else
    cli::cli_alert_warning("Copy Number data has not been used for this analysis.")

  n = sum(x$data$OK_tumour)
  N = nrow(x$data)
  
  cli::cli_alert(
    "Mutations data: n = {.field {n}} out of {.value {N}} within range ({.field {round(n/N, 2) * 100}%})."
  )
  cat("\n")

  pio::pioStr('        TIT : ',
              paste0(round(x$TIT, 2) * 100, '% (RF ', round(x$TIT_rf, 2) * 100, ')'),
              ' ~ n =', length(x$fit$mobster_analysis$clonal_mutations),
              'clonal mutations, cluster', x$fit$mobster_analysis$clonal_cluster,
              suffix = '\n')

  pio::pioStr('        TIN : ',
              paste0(round(x$TIN, 2) * 100, '% (RF ', round(x$TIN_rf, 2) * 100, ')'),
              ' ~ n =',
              sum(x$fit$BMix_analysis$output$VAF > 0), 'with VAF > 0',
              suffix = '\n')

  # Report Classification
  cat('\n')

  cn = classification_normal(x$TIN)
  ct = classification_tumour(x$TIT)

  # pc = c('forestgreen', 'steelblue', 'goldenrod1', 'indianred3', 'purple')

  if(ct['color'] == 'forestgreen') cat(black$bgGreen$bold("   QC Tumour  "), ct['QC'] %>% green)
  if(ct['color'] == 'steelblue') cat(black$bgBlue$bold("   QC Tumour  "), ct['QC'] %>% blue)
  if(ct['color'] == 'goldenrod4') cat(black$bgYellow$bold("   QC Tumour  "), ct['QC'] %>% yellow)
  if(ct['color'] == 'indianred3') cat(black$bgRed$bold("   QC Tumour  "), ct['QC'] %>% red)
  if(ct['color'] == 'purple') cat(black$bgMagenta$bold("   QC Tumour  "), ct['QC'] %>% magenta)

  cat('\n')

  if(cn['color'] == 'forestgreen') cat(black$bgGreen$bold("   QC Normal  "), cn['QC'] %>% green)
  if(cn['color'] == 'steelblue') cat(black$bgBlue$bold("   QC Normal  "), cn['QC'] %>% blue)
  if(cn['color'] == 'goldenrod4') cat(black$bgYellow$bold("   QC Normal  "), cn['QC'] %>% yellow)
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
#' @description This function is a wrapper to call \code{plot_simple_report}.
#'
#' @seealso \code{plot_simple_report}.
#'
#' @param x A TINC analysis computed with \code{autofit}.
#' @param ... Extra S3 parameters.
#'
#' @return A TINC analysis plot (a `ggplot` figure with multiple panels).
#' @export
#'
#' @examples
#' rt = random_TIN()
#' plot(autofit(input = rt$data, cna = rt$cna, FAST = TRUE))
plot.tin_obj = function(x, ...)
{
  if(!inherits(x, "tin_obj")) stop("Not a TINC object .... run autofit(.) first, aborting.")

  # plot_full_page_report(x$fit)
  plot_simple_report(x)
}

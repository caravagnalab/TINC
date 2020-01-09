#' Title
#'
#' @param N
#' @param TIN
#' @param TIT
#' @param normal_coverage
#' @param tumour_coverage
#'
#' @return
#' @export
#'
#' @examples
#' random_TIN()
random_TIN = function(N = 1000, TIN = 0.05, TIT = 1, normal_coverage = 30, tumour_coverage = 120)
{
  # SNVs
  mb = mobster::random_dataset(N, K_betas = 1, Beta_bounds = c((TIT/2)-0.01, (TIT/2) + 0.01),
                               Beta_variance_scaling = runif(1, 200, 500))

  factor_correction = ((TIN/2)/(TIT/2))

  L_sample = pio:::nmfy(
    CNAqc::chr_coordinates_hg19$chr,
    CNAqc::chr_coordinates_hg19$length
  )

  z = mb$data
  z$chr = paste0('chr', sample(1:22, N, replace = T))
  z$from = sapply(z$chr, function(x) sample(1:L_sample[x], 1))
  z$to = z$from + 1
  z$ref = sample(c("A", "C", "T", "G"), replace = TRUE, nrow(z))
  z$alt = sample(c("A", "C", "T", "G"), replace = TRUE, nrow(z))

  z = z %>% dplyr::select(chr, from, to, ref, alt, VAF, cluster) %>%
    dplyr::mutate(
      T.VAF = VAF,
      N.VAF = T.VAF * factor_correction
    )

  z$n_tot_count = rpois(nrow(z), lambda = normal_coverage)
  z$t_tot_count = rpois(nrow(z), lambda = tumour_coverage)

  z$n_alt_count = floor(z$n_tot_count * z$N.VAF)
  z$t_alt_count = floor(z$t_tot_count * z$T.VAF)

  z$n_ref_count = z$n_tot_count - z$n_alt_count
  z$t_ref_count =  z$t_tot_count - z$t_alt_count

  z$sim_t_vaf =  z$t_alt_count/z$t_tot_count
  z$sim_n_vaf =  z$n_alt_count/z$n_tot_count

  cli::cli_alert_success(
    'Generated TINC dataset (n = {.value {N}} mutations), TIN ({.value {TIN}}) and TIT ({.value {TIT}}), normal and tumour coverage {.value {normal_coverage}}x and {.value {tumour_coverage}}x.'
    )

  plot =
    ggplot(z, aes(x = sim_n_vaf, y = sim_t_vaf)) +geom_point()  +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(title = "TIN simulated data",
         x = 'Simulated normal VAF',
         y = 'Simulated tumour VAF',
         subtitle = paste0('TIN = ', TIN, ' ~ TIT = ', TIT)) +
    geom_vline(xintercept = TIN/2, linetype = 'dashed', size = .3) +
    geom_hline(yintercept = TIT/2, linetype = 'dashed', size = .3) +
    mobster:::my_ggplot_theme()

  # CNA
  td_cna =   TINC:::as_tumour(z) %>%
    dplyr::mutate(minor = 1, Major = 1, from = from - 1, to = to + 1) %>%
    dplyr::select(-ref, -alt, -DP, -NV, -VAF)

  return(
    list(
      data =   z %>%
        dplyr::select(-VAF, -T.VAF, -N.VAF),
      cna = td_cna,
      plot = plot
    )
  )
}

random_CNA = function(x)
{
  CNAqc::chr_coordinates_hg19 %>%
    dplyr::select(-length, -starts_with('centromer')) %>%
    dplyr::mutate(minor = 1, Major = 1)
}

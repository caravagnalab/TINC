#' Title
#'
#' @param N
#'
#' @return
#' @export
#'
#' @examples
#' random_TIN()
random_TIN = function(N = 1000, TIN = 0.05, TIT = 1)
{
  mb = mobster::random_dataset(N, K_betas = 1, Beta_bounds = c((TIT/2)-0.01, (TIT/2) + 0.01),
                               Beta_variance_scaling = runif(1, 200, 500))

  factor_correction = ((TIN/2)/(TIT/2))

  z = mb$data %>%
    mutate(
      id = paste('chr_fake', row_number(), row_number() + 1, 'A', 'C', sep = ':'),
      filters = 'PASS',
      T.VAF = VAF,
      N.VAF = T.VAF * factor_correction
    )

  z$n_tot_count = rpois(nrow(z), lambda = 30)
  z$t_tot_count = rpois(nrow(z), lambda = 120)

  z$n_alt_count = floor(z$n_tot_count * z$N.VAF)
  z$t_alt_count = floor(z$t_tot_count * z$T.VAF)

  z$n_ref_count = z$n_tot_count - z$n_alt_count
  z$t_ref_count =  z$t_tot_count - z$t_alt_count

  z$sim_t_vaf =  z$t_alt_count/z$t_tot_count
  z$sim_n_vaf =  z$n_alt_count/z$n_tot_count

  # ggplot(z, aes(x = N.VAF, y = sim_n_vaf)) + geom_point() + geom_abline()
  # ggplot(z, aes(x = T.VAF, y = sim_t_vaf)) + geom_point() + geom_abline()

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

  z$filters = "PASS"
  return(
    list(
      data =   z %>%
        dplyr::select(id, n_ref_count, n_alt_count, t_ref_count, t_alt_count, filters),
      plot = plot
    )
  )
}


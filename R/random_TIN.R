#' Generate a random dataset for TINC.
#' 
#' @description This function simulates a random dataset for TINC analysis,
#' with mutations and copy number data. Segments are not real, they are assumed
#' to be constant heterozygous diploid (`Major = minor = 1`) and span just
#' mutations (for mappability).
#' 
#' This samples has some noise so that the obtained TIT score might be slightly
#' lower than the required input.
#'
#' @param N Number of input simulations
#' @param TIN TIN - Tumour in normal contamination level.
#' @param TIT TIT - Tumour in tumour contamination level (aka tumour purity).
#' @param normal_coverage Normal coverage (mean). 
#' @param tumour_coverage Tumour coverage (mean).
#'
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_hline geom_vline labs xlim ylim theme_light theme geom_point
#' @return Tibbles with the data, and a plot.
#' @export
#' @md
#'
#' @examples
#' 
#' set.seed(1234)
#' 
#' # Default dataset
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

  z = mb$data %>% filter(VAF < ((TIT/2) * 1.5))
  
  N = z %>% nrow
  
  z$chr = paste0('chr', sample(1:22, N, replace = T))
  z$from = sapply(z$chr, function(x) sample(1:L_sample[x], 1))
  z$to = z$from + 1
  z$ref = sample(c("A", "C", "T", "G"), replace = TRUE, nrow(z))
  z$alt = sample(c("A", "C", "T", "G"), replace = TRUE, nrow(z))

  z = z %>% dplyr::select(chr, from, to, ref, alt, VAF, simulated_cluster) %>%
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
  
  p_t = z %>% 
    ggplot() + 
    geom_histogram(aes(T.VAF), binwidth = 0.01) +
    my_ggplot_theme() +
    xlim(0,1) +
    labs(x = "Tumour VAF")
  
  p_n =   z %>% 
    ggplot() + 
    geom_histogram(aes(N.VAF), binwidth = 0.01) +
    my_ggplot_theme() +
    xlim(0,1) +
    labs(x = "Normal VAF")

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
    my_ggplot_theme()
  
  plot = 
    cowplot::plot_grid(
      plot,
      cowplot::plot_grid(p_t, p_n, nrow = 2),
      axis = 'tb', 
      align = 'h'
    )

  # CNA
  td_cna =  as_tumour(z) %>%
    dplyr::mutate(minor = 1, Major = 1, from = from - 1, to = to + 1) %>%
    dplyr::select(chr, from, to, ref, Major, minor)

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
    dplyr::mutate(minor = 1, Major = 1) %>%
    dplyr::mutate(length = to - from)
}

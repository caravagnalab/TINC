#' Runs a TINC analysis.
#'
#' @description This function is a wrapper to run the main analysis of TINC.
#'
#' The steps are as follows:
#' \itemize{
#' \item 1) Input data is loaded from a `file`, or from a `dataframe`.
#' \item 2) Clonal mutations are estimated for the tumour, together with
#' the tumour purity (Tumour in Tumour).
#' \item 3) From putative clonal mutations of the tumour, the Tumour in Normal
#' contamination level is estimated.
#' }
#' An S3 object is returned that contains the results of the analysis.
#'
#' @param input A `file` to load from disk, or a `dataframe `consistent
#' with the required input. See the vignette for more information.
#' @param VAF_range_tumour A range `[x, y]` so that only mutations
#' with VAF in that range are actually used to determine the TIN/ TIT
#' levels of the input.
#' @param cutoff_miscalled_clonal An upper bound on the VAF of a
#' cluster in the tumour data. Clusters above this value will be
#' considered miscalled clonal clusters (e.g., due to LOH etc.).
#' @param cutoff_lv_assignment Consider only latent variables with
#' responsibilities above this cutoff.
#' @param N If there are more than `N` mutations in VAF range
#' `VAF_range_tumour`, a random subset of size `N` is retained.
#' @param FAST If `TRUE`, it runs the analysis with reduced sampling
#' power and accuracy. Use this to obtain a result for preliminary
#' inspection of your data, and then run `autofit` with this parameter
#' set to `FALSE`.
#'
#' @return An S3 object that contains the results of this analysis.
#'
#' @export
#'
#' @import dplyr
#' @import mobster
#' @import BMix
#' @import VIBER
#'
#' @examples
#' autofit(random_TIN(), FAST = TRUE)
autofit = function(input,
                   VAF_range_tumour = c(0, 0.7),
                   cutoff_miscalled_clonal = .6,
                   cutoff_lv_assignment = 0.75,
                   N = 20000,
                   FAST = FALSE)
{
  # Load data, requires the "VAF_range_tumour" status
  dataset = load_input(input, VAF_range_tumour = VAF_range_tumour, N = N)

  # Plot input data
  # dataset_plot_raw = plot_raw(dataset = dataset)

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # MOBSTER fit of the tumour. It fits the tumour, determines the clonal cluster,
  # a pool of highly-confidence clonal mutations and estimates the purity of the tumour
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # mobster_analysis = wrap_mobster(tumour, cutoff_miscalled_clonal, cutoff_lv_assignment, ...)

  pio::pioTit("MOBSTER fit", 'n =', nrow(dataset$tumour %>% filter(used)))

  if(FAST)
    dataset$mobster_analysis = analyse_mobster(
      dataset$tumour %>% filter(used),
      cutoff_miscalled_clonal,
      cutoff_lv_assignment,
      auto_setup = 'FAST'
    )
  else
    dataset$mobster_analysis = analyse_mobster(
      dataset$tumour %>% filter(used),
      cutoff_miscalled_clonal,
      cutoff_lv_assignment,
      K = 1:3,
      samples = 6,
      maxIter = 300,
      parallel = FALSE,
      epsilon = 1e-9,
      init = 'random'
    )

  # pio::pioStr(
  #   "\n  Clonal cluster",
  #   dataset$mobster_analysis$clonal_cluster,
  #   '~ n =',
  #   mobster_fit_tumour$best$N.k[clonal_cluster],
  #   suffix = '\n'
  # )
  # pio::pioStr("Estimated purity", estimated_tumour_purity, suffix = '\n')

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # BMix fit of the germline read counts of clonal mutations identified by MOBSTER.
  # It fits the tumour, determines the clonal cluster,
  # a pool of highly-confidence clonal mutations and estimates the purity of the tumour
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  g_input = dataset$normal %>%
    dplyr::filter(id %in% dataset$mobster_analysis$clonal_mutations)

  pio::pioTit("BMix fit", 'n =', nrow(g_input), prefix = '\t')

  # Normal sample ~ use putative clonal mutations from MOBSTER
  if(FAST)
  dataset$BMix_analysis = analyse_BMix(
    x = g_input,
    K.BetaBinomials = 0,
    epsilon = 1e-6,
    samples = 2
    )
  else
    dataset$BMix_analysis = analyse_BMix(
      x = g_input,
      K.BetaBinomials = 0,
      epsilon = 1e-9,
      samples = 6
    )

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # TIN profiling (plus TIT)
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  final_TIN_profile = TIN_profiler(dataset)

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # VIBER profiling
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pio::pioTit("VIBER fit", 'n =', nrow(dataset$joint), prefix = '\t')

  if(FAST)
    final_TIN_profile$VIBER_analysis = analyze_VIBER(dataset$joint,
                                                     K = 5,
                                                     alpha_0 = 1e-6,
                                                     max_iter = 1000,
                                                     epsilon_conv = 1e-6,
                                                     samples = 3)
  else
    final_TIN_profile$VIBER_analysis = analyze_VIBER(dataset$joint,
                                                     K = 8,
                                                     alpha_0 = 1e-6,
                                                     max_iter = 5000,
                                                     epsilon_conv = 1e-9,
                                                     samples = 10)

  # Output TINC object of class tin_obj
  output_obj = list(
    fit = final_TIN_profile,
    TIN = final_TIN_profile$BMix_analysis$estimated_purity,
    TIT = final_TIN_profile$mobster_analysis$estimated_purity,
    params = list(
      VAF_range_tumour = VAF_range_tumour,
      cutoff_miscalled_clonal = cutoff_miscalled_clonal,
      cutoff_lv_assignment = cutoff_lv_assignment,
      N = N,
      FAST = FAST
    )
  )

  class(output_obj) = "tin_obj"

  return(output_obj)
}

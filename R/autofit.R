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
                   cutoff_lv_assignment = 0,
                   N = 20000,
                   FAST = FALSE)
{
  # Load data, requires the "VAF_range_tumour" status
  dataset = load_input(input, VAF_range_tumour = VAF_range_tumour, N = N)

  # Plot input data
  dataset_plot_raw = plot_raw(dataset = dataset)

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
  dataset$BMix_analysis = analyse_BMix(
    x = g_input,
    K.BetaBinomials = 0)

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # TIN profiling (plus TIT)
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  final_TIN_profile = TIN_profiler(dataset)

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # VIBER profiling
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  pio::pioTit("VIBER fit", 'n =', nrow(dataset$joint), prefix = '\t')
  final_TIN_profile$VIBER_analysis = analyze_VIBER(dataset$joint)

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Report assembly
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Plots of all the data
  data_fit_panel = cowplot::plot_grid(
    dataset_plot_raw,
    plot_sample_contamination(final_TIN_profile),
    rel_widths = c(1, .7),
    align = 'h',
    # axis = 'bt',
    labels = c("A", "B"),
    nrow = 1,
    ncol = 2
  )

  # data_fit_panel =
  #   ggpubr::annotate_figure(data_fit_panel,
  #                           top = ggpubr::text_grob(
  #                             paste0('Input: n = ', nrow(file, '\n'),
  #                             hjust = 0,
  #                             x = 0,
  #                             size = 18
  #                           ))

  fit_panel =
    cowplot::plot_grid(
      final_TIN_profile$mobster_analysis$plot,
      final_TIN_profile$BMix_analysis$plot,
      nrow = 1,
      ncol = 2,
      align = 'h',
      axis = 'bt',
      labels = c('C', 'D')
    )

  sq_panel = cowplot::plot_grid(
    plot_contamination_full_size(final_TIN_profile),
    plot_contamination_zoom(final_TIN_profile),
    final_TIN_profile$VIBER_analysis$plot,
    nrow = 1,
    ncol = 3,
    align = 'h',
    axis = 'bt',
    rel_widths = c(2, 1, 1),
    labels = c('E', '', 'D')
  )

  full_figure =  cowplot::plot_grid(
    data_fit_panel,
    fit_panel,
    sq_panel,
    nrow = 3,
    ncol = 1,
    align = 'v',
    axis = 'lr',
    rel_heights = c(.8, .7, .7)
  )

  output_obj = list(
    fit = final_TIN_profile,
    plot = full_figure,
    TIN = final_TIN_profile$BMix_analysis$estimated_purity,
    TIT = final_TIN_profile$mobster_analysis$estimated_purity
  )

  class(output_obj) = "tin_obj"

  return(output_obj)
}

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
#' @param input A`dataframe` of the iput mutations. Must be in a certain format,
#' see the vignette for more information.
#' @param cna Copy Number data in the format of package \code{CNAqc}.
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
#' # Random
#' rt = TINC::random_TIN()
#' x = TINC::autofit(input = rt$data, cna = rt$cna, FAST = TRUE)
#' print(x)
autofit = function(input,
                   cna,
                   VAF_range_tumour = c(0, 0.7),
                   cutoff_miscalled_clonal = .6,
                   cutoff_lv_assignment = 0.75,
                   N = 20000,
                   FAST = FALSE)
{
  mobster_analysis = BMix_analysis = VIBER_analysis = NULL

  pio::pioHdr("TINC")
  cat('\n')

  # Load data, checks VAF range and N limit; subset with CNA data if required
  input_data = TINC:::load_TINC_input(input, cna, VAF_range_tumour = VAF_range_tumour, N = N)

  x = input_data$mutations
  cna_map = input_data$cna_map

  # Used karyptype, only for CNA-based analyses
  used_karyotype = ifelse(
    TINC:::analysis_mode(cna) != "CNA",
    NA,
    input_data$what_we_used
  )

  # By default, we assume we are looking at the whole genome. When CNA are available, we restrict the set
  used_chromosomes = paste0('chr', 1:22)
  if(TINC:::analysis_mode(cna) == "CNA")
    used_chromosomes = unique(cna$chr)

  # TODO maybe used_chromosomes should be subset according to used_karyotype?

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # MOBSTER fit of the tumour. It fits the tumour, determines the clonal cluster,
  # a pool of highly-confidence clonal mutations and estimates the purity of the tumour
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  if (FAST)
    mobster_analysis = TINC:::analyse_mobster(
      x = TINC:::as_tumour(x) %>% dplyr::filter(OK_tumour),
      cna_map = cna_map,
      cutoff_miscalled_clonal = cutoff_miscalled_clonal,
      cutoff_lv_assignment = cutoff_lv_assignment,
      chromosomes = used_chromosomes,
      auto_setup = 'FAST'
    )
  else{
    # When we are asked not to run the fast analysis, we actually decide what is a reasonable
    # value of K. If we have CNA data, we use that as a strong "prior" (not in a Bayesian sense);
    # otherwise we just use the default
    K = 1:3
    if(TINC:::analysis_mode(cna) == "CNA"){
      # The prior is on the mixture expected
      if(used_karyotype %in% c("2:0", "2:1", "2:2")) K = 2:3
    }

    mobster_analysis = TINC:::analyse_mobster(
      x = TINC:::as_tumour(x) %>% dplyr::filter(OK_tumour),
      cna_map = cna_map,
      cutoff_miscalled_clonal = cutoff_miscalled_clonal,
      cutoff_lv_assignment = cutoff_lv_assignment,
      chromosomes = used_chromosomes,
      K = K,
      samples = 6,
      maxIter = 300,
      parallel = FALSE,
      epsilon = 1e-9,
      init = 'random'
    )
  }

  mobster_analysis_dynamic_cutoff = (cutoff_lv_assignment != mobster_analysis$cutoff_lv_assignment)

  cat('\n')
  x$OK_clonal = x$id %in% mobster_analysis$clonal_mutations
  cli::cli_alert_success(
    "MOBSTER found n = {.value {sum(x$OK_clonal)}} clonal mutations from cluster {.value {mobster_analysis$clonal_cluster}}"
  )

  # cli::cli_process_done()


  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # BMix fit of the germline read counts of clonal mutations identified by MOBSTER.
  # It fits the tumour, determines the clonal cluster,
  # a pool of highly-confidence clonal mutations and estimates the purity of the tumour
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  # cli::cli_process_start("BMix analysis of normal data")

  # Normal sample ~ use putative clonal mutations from MOBSTER
  if (FAST)
    BMix_analysis = TINC:::analyse_BMix(
      x = TINC:::as_normal(x) %>%
        dplyr::filter(OK_clonal),
      cna_map = cna_map,
      K.BetaBinomials = 0,
      epsilon = 1e-6,
      samples = 2
    )
  else
    BMix_analysis = TINC:::analyse_BMix(
      x = as_normal(x) %>%
        dplyr::filter(OK_clonal),
      cna_map = cna_map,
      K.BetaBinomials = 0,
      epsilon = 1e-9,
      samples = 6
    )

  # cli::cli_process_done()

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # VIBER profiling
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # pio::pioTit("VIBER fit", 'n =', nrow(dataset$joint), prefix = '\t')

  # cli::cli_process_start("VIBER analysis of joint data")

  if (FAST)
    VIBER_analysis = TINC:::analyze_VIBER(
      x = x,
      K = 5,
      alpha_0 = 1e-6,
      max_iter = 1000,
      epsilon_conv = 1e-6,
      samples = 3
    )
  else
    VIBER_analysis = analyze_VIBER(
      x = x,
      K = 8,
      alpha_0 = 1e-6,
      max_iter = 5000,
      epsilon_conv = 1e-9,
      samples = 10
    )

  # cli::cli_process_done()

  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  # Full profiling
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  x = x %>%
    full_join(mobster_analysis$output %>% dplyr::select(id, tumour.cluster),
              by = c("id")) %>%
    full_join(BMix_analysis$output %>% dplyr::select(id, normal.cluster),
              by = c("id")) %>%
    full_join(VIBER_analysis$output,
              by = c("id"))


  # Output TINC object of class tin_obj
  cna_list = list(NULL)
  if(analysis_mode(cna) == "CNA")
    {
      cna_list = data.frame(
        used_chromosomes = paste(used_chromosomes, collapse = ':'),
        karyotype = input_data$what_we_used,
        ploidy = cna_map$ploidy,
        stringsAsFactors = FALSE
      )
  }




  # Output object
  output_obj = list(
    data = x,
    analysis_type = TINC:::analysis_mode(cna),
    fit = list(
      BMix_analysis = BMix_analysis,
      mobster_analysis = mobster_analysis,
      VIBER_analysis = VIBER_analysis,
      CNA = cna_map
    ),
    TIN = BMix_analysis$estimated_purity,
    TIT = mobster_analysis$estimated_purity,
    TIN_rf = BMix_analysis$estimated_read_fraction,
    TIT_rf = mobster_analysis$estimated_read_fraction,
    params = list(
      VAF_range_tumour = VAF_range_tumour,
      cutoff_miscalled_clonal = cutoff_miscalled_clonal,
      cutoff_lv_assignment = cutoff_lv_assignment,
      cutoff_lv_assignment_dynamic = mobster_analysis_dynamic_cutoff,
      N = N,
      FAST = FAST
    ),
    cna_params = cna_list
  )

  class(output_obj) = "tin_obj"

  return(output_obj)
}

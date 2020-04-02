#' Load TINC input data.
#'
#' @description
#'
#' The function loads a dataframe or tibble, with a set of required
#' column names. See the package for a detailed description of the
#' input format.
#'
#' After loading data, mutation with VAF outside a range (default [0; 0.7])
#' are flagged and removed from downnstream analysis.
#'
#' @param x A dataframe or tibble.
#' @param cna Copy Number data in the format of package \code{CNAqc}.
#' @param VAF_range_tumour 2D vector for a VAF range used to filter
#' mutations from the tumour sample.
#' @param N If there are more than `N` mutations in VAF range
#' `VAF_range_tumour`, a random subset of size `N` is retained.
#'
#' @return A tibble with the loaded data.
#'
#' @examples
#' # Generating a random TIN input
#' load_TINC_input(random_TIN())
load_TINC_input = function(x,
                           cna,
                           VAF_range_tumour = c(0, 0.7),
                           N = 20000)
{
  cli::cli_h1("Loading TINC input data")

  stopifnot(is.data.frame(x))

  # Check columns
  required_colnames = c(
    'chr',
    'from',
    'to',
    'ref',
    'alt',
    'n_ref_count',
    'n_alt_count',
    't_ref_count',
    't_alt_count'
  )

  stopifnot(all(required_colnames %in% colnames(x)))

  # Data
  x = x[, required_colnames] %>%
    as_tibble() %>%
    dplyr::mutate(
      karyotype = "Unknown",
      id = paste(chr, from, to, ref, alt, sep = ':')
    )

  most_prevalent_karyotype = NULL

  # cli::cli_alert_info("Using for mutation data {.url {required_colnames}} for n = {.value {nrow(x)}}.")
  cli::cli_alert_success("Input data contains n = {.value {nrow(x)}} mutations, selecting operation mode.")

  ####################################
  # With CNA data: establish a special execution setup of this run
  ####################################
  cn_obj = what_we_used = NULL
  if(analysis_mode(cna) == "CNA")
  {
    cli::cli_alert_warning("Found CNA data, retaining only mutations that map to segments with predominant karyotype ...")

    # Map CNA data, and retain only mappable mutations that belong to the largest karyotype
    cat("\n")
    cn_obj = CNAqc::init(snvs = TINC:::as_tumour(x), cna = cna, .8)
    cat("\n")

    # CNAqc computes the most prevalent karyotype with the largest prevalence (in basepairs).
    most_prevalent_karyotype = cn_obj$most_prevalent_karyotype

    supported_karyotypes = c('1:0', '1:1', '2:1', '2:2', '2:0')

    cli::cli_h2("Genome coverage by karyotype, in basepairs.\n")
    print(cn_obj$basepairs_by_karyotype)

    if(!(most_prevalent_karyotype %in% supported_karyotypes))
    {
      cli::cli_alert_danger(
        "The most prevalent karyotype (in basepairs) is {.field {most_prevalent_karyotype}} and is not any of {.field {supported_karyotypes}}, will use the largest among those instead ..."
        )

      most_prevalent_karyotype = cn_obj$basepairs_by_karyotype %>%
        dplyr::filter(karyotype %in% supported_karyotypes) %>%
        dplyr::filter(row_number() == 1) %>%
        dplyr::pull(karyotype)

      if(length(most_prevalent_karyotype) == 0) {
        stop("We did not find any of the supported karyotypes, will not analyse this samples with CNA data.")
      }
    }

    # Extract those mutations
    mappable = cn_obj$snvs %>% dplyr::filter(karyotype == most_prevalent_karyotype) %>% pull(id)

    x = x %>%
      dplyr::filter(id %in% mappable) %>%
      dplyr::mutate(karyotype = most_prevalent_karyotype)

    cli::cli_alert_success("n = {.value {nrow(x)}} mutations mapped to CNA segments with karyotype {.field {most_prevalent_karyotype}} (largest available in basepairs).")
  }

  # Tumour filter
  tumour_exclude_VAF_range = as_tumour(x) %>%
    dplyr::filter((VAF < VAF_range_tumour[1]) |
                    (VAF > VAF_range_tumour[2])) %>%
    dplyr::pull(id)

  x = x %>%
    dplyr::mutate(# OK_tumour = TRUE,
      OK_tumour = !(id %in% tumour_exclude_VAF_range))

  cli::cli_alert_success(
    "Mutation with VAF within {.value {VAF_range_tumour}} ~ n = {.value { sum(x$OK_tumour)}}."
  )

  # Downsample
  if (sum(x$OK_tumour) > N)
  {
    cli::cli_alert_warning("More than {.value {N}} mutations, downsampling.")

    w_t =  x %>% dplyr::filter(OK_tumour) %>% dplyr::sample_n(N) %>% dplyr::pull(id)
    x$OK_tumour[!(x$id %in% w_t)] = FALSE
  }

  return(list(mutations = x, cna_map = cn_obj, what_we_used = most_prevalent_karyotype))
}


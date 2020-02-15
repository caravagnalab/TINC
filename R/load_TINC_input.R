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

  # cli::cli_alert_info("Using for mutation data {.url {required_colnames}} for n = {.value {nrow(x)}}.")
  cli::cli_alert_success("Input data contains n = {.value {nrow(x)}} mutations, selecting operation mode.")

  ####################################
  # With CNA data: establish a special execution setup of this run
  ####################################
  cn_obj = NULL
  if(analysis_mode(cna) == "CNA")
  {
    cli::cli_alert_warning("Found CNA data, retaining only mutations that map to segments with predominant karyotype ...")

    # Map CNA data, and retain only mappable mutations that belong to the largest karyotype
    cat("\n")
    cn_obj = CNAqc::init(snvs = TINC:::as_tumour(x), cna = cna, .8)
    cat("\n")

    # First mappable in general
    mappable = cn_obj$snvs %>% dplyr::filter(!is.na(segment_id))

    # CNAqc estimates the ploidy as the karyotypes with the largest prevalence (in basepairs).
    #
    # We subset CNAs according to the prevalent karyotype (defined from total ploidy, p). For p = 1 we take LOH regions,
    # for p = 2 we take balanced diploid regions (1:1), for p = 3 triploid (2:1), for p = 4 tetraploid (2:2).
    # For anything else, we take diploid ones.
    ploidy_by_karyotypes = round(cn_obj$ploidy)

    if(ploidy_by_karyotypes == 1) mappable = mappable %>% dplyr::filter(karyotype =='1:0')
    if(ploidy_by_karyotypes == 2) mappable = mappable %>% dplyr::filter(karyotype =='1:1')
    if(ploidy_by_karyotypes == 3) mappable = mappable %>% dplyr::filter(karyotype =='2:1')
    if(ploidy_by_karyotypes == 4) mappable = mappable %>% dplyr::filter(karyotype =='2:2')

    if(ploidy_by_karyotypes > 4) {

      largest_karyotypes_below_4 = cn_obj$n_karyotype[which(names(cn_obj$n_karyotype) %in% c('1:0', '1:1', '2:1', '2:2'))] %>% sort(decreasing = TRUE)
      largest = largest_karyotypes_below_4[1] %>% names

      cli::cli_alert_danger("The ploidy of this tumour is > 4, we will use as karyotype {.field {largest}} which is the largest one with ploidy <= 4.")

      mappable = mappable %>% dplyr::filter(karyotype == largest)
    }

    what_we_used = mappable$karyotype[1]

    cli::cli_alert_success("n = {.value {nrow(x)}} mutations mapped to {.value {nrow(cna)}} CNA segments with karyotype {.field {what_we_used}} (largest available).")

    mappable = mappable %>% dplyr::pull(id)
    x = x %>%
      dplyr::filter(id %in% mappable) %>%
      dplyr::mutate(karyotype = what_we_used)
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

  return(list(mutations = x, cna_map = cn_obj))
}


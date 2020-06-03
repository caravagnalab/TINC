#' Canvas VCF parsing function
#'
#' @description Parse a VCF file from Canvas, which stores a
#' list of segments with absolute CNA events, and tumour purity.
#'
#' The obtained data, together with mutation calls (loaded elsewhere)
#' can be used to call function \code{autofit}.
#'
#' @param file VCF filename.
#'
#' @return A named list with the calls and tumour purity.
#'
#' @export
#'
#' @importFrom vcfR read.vcfR
#'
#' @examples
#' # not run
#'
#'
load_VCF_Canvas = function(file) {
  if (!file.exists(file))
    stop("Input file", file, "not found!")

  cli::cli_h2("VCF Canvas parser for TINC")

  # Read VCF with vcfR
  segments = vcfR::read.vcfR(file, verbose = FALSE)

  # Genotypes
  gt = vcfR::extract_gt_tidy(segments) %>%
    dplyr::mutate(
      Major = gt_MCC,
      minor = gt_CN - Major
    ) %>%
    select(minor, Major)

  # Locations
  locs = segments@fix[, "ID"] %>%
    dplyr::as_tibble() %>%
    tidyr::separate(value, into = c("A", "SD", "chr", "range"), sep = ":") %>%
    tidyr::separate(range, into = c("from", "to"), sep = "-") %>%
    dplyr::mutate(
      from = as.numeric(from),
      to = as.numeric(to)
      ) %>%
    select(-A, -SD)

  # Canvas calls
  canvas_calls = dplyr::bind_cols(locs, gt)
  canvas_calls = canvas_calls[complete.cases(canvas_calls), ] %>%
    mutate(length = to - from)

  # Purity is in the VCF as well
  tumour_purity = strsplit(
    unlist(
      vcfR::queryMETA(segments, element = 'EstimatedTumorPurity', nice = TRUE)
    ),
    split = "="
  )[[1]][2]

  tumour_purity = as.numeric(tumour_purity)

  # Report some text
  cli::cli_alert_success(
    "Canvas calls parsed successfully: n = {.value {nrow(canvas_calls)}} CNA segments."
  )

  print(canvas_calls)

  cli::cli_alert_success(
    "Canvas tumour purity: p = {.value {tumour_purity}}."
  )

  return(list(calls = canvas_calls, purity = tumour_purity))
}

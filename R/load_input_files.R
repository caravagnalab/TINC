#' Load data for TIN.
#'
#' @description
#'
#' The function can support both a file argument and a dataframe.
#' In the case a file is given, no column or row names should be
#' used. If the input is a dataframe or tibble, column names are
#' required. See the package for a detailed description of the
#' input format for TIN analysis.
#'
#' This function, after loading data, annotates which mutation
#' has VAF  within a custom range (default [0; 0.7]); these are
#' the only mutations that are used for the TIN analysis.
#'
#' @param x Either a filename, or a dataframe.
#' @param verbose If `TRUE`, more output is print to screen.
#' @param VAF_range_tumour A 2D vector that determines a
#' VAF range used to decide which mutation we do analyse from the
#' input cancer sample.
#' @param N If there are more than `N` mutations in VAF range
#' `VAF_range_tumour`, a random subset of size `N` is retained.
#'
#' @return A list that contains the following informations:
#' tumour calls, normal calls and joint calls, plus the input
#' VAF range and file name.
#'
#' @export
#'
#' @examples
#' # Generating a random TIN input
#' load_input(random_TIN())
load_input = function(x,
                      verbose = FALSE,
                      VAF_range_tumour = c(0, 0.7),
                      N = 20000)
{
  file = x
  pio::pioHdr("TINC ~ Load tumour/normal data")

  ######## Detect if it is filename or data.frame, and open it
  if (is.character(file) && file.exists(file))
  {
    pio::pioStr(
      "\n    Input file ",
      file,
      '~',
      file.size(file) %>% utils:::format.object_size(units = "auto"),
      suffix = '\n'
    )

    # VCF file dumped in our TIN format
    dataset =  read.csv(
      file,
      sep = '\t',
      header = FALSE,
      stringsAsFactors = FALSE
    ) %>% as_tibble

    if (ncol(dataset) > 5) {
      message("Input VCF file has more than 5 columns, only the first 5 will be used.")
      dataset = dataset[, 1:5]
    }

    colnames(dataset) = c('id',
                          'n_ref_count',
                          'n_alt_count',
                          't_ref_count',
                          't_alt_count')

  }

  if (is.data.frame(file))
  {
    pio::pioStr("\n    Input dataframe ", dim(file)[1],  'x', dim(file)[2],
                suffix = '\n')

    # required TIN format
    required_cols = c('id',
                      'n_ref_count',
                      'n_alt_count',
                      't_ref_count',
                      't_alt_count')

    if (!all(required_cols %in% colnames(file)))
      stop(
        "The input dataframe must have the following named columns: ",
        paste(required_cols, collapse = ', ')
      )

    dataset = file %>% as_tibble()
  }
  ######## Detected: DONE


  # Germline data
  germline = dataset %>%
    dplyr::select(id, starts_with('n')) %>%
    dplyr::mutate(DP = n_ref_count + n_alt_count,
                  NV = n_alt_count,
                  VAF = NV / DP) %>%
    tidyr::separate(
      id,
      into = c('chr', 'from', 'to', 'ref', 'alt'),
      sep = ':',
      remove = FALSE
    ) %>%
    dplyr::select(chr,
                  from,
                  to,
                  ref,
                  alt,
                  id,
                  ends_with('count'),
                  DP,
                  NV,
                  VAF)

  pio::pioTit("Mutations annotated in the normal sample")
  pio::pioDisp(germline)

  #  Tumour data
  tumour = dataset %>%
    dplyr::select(id, starts_with('t')) %>%
    dplyr::mutate(DP = t_ref_count + t_alt_count,
                  NV = t_alt_count,
                  VAF = NV / DP) %>%
    tidyr::separate(
      id,
      into = c('chr', 'from', 'to', 'ref', 'alt'),
      sep = ':',
      remove = FALSE
    ) %>%
    dplyr::select(chr,
                  from,
                  to,
                  ref,
                  alt,
                  id,
                  ends_with('count'),
                  DP,
                  NV,
                  VAF)

  pio::pioTit("Mutations annotated in the tumour sample")
  pio::pioDisp(tumour)

  pio::pioStr("\n VAF summaries ", suffix = '\n\n')
  bind_rows(
    c(`sample` = "Normal", summary(germline$VAF) %>% round(2)),
    c(`sample` = "Tumour", summary(tumour$VAF) %>% round(2))
  ) %>%
    print

  # Joint data
  joint_data = germline %>%
    dplyr::full_join(
      tumour,
      by = c('chr', 'from', 'to', 'ref', 'alt', 'id'),
      suffix = c('.normal', '.tumour')
    )

  # Tumour filter
  tumour = tumour %>%
    dplyr::mutate(used = (VAF > VAF_range_tumour[1]) &
                    (VAF < VAF_range_tumour[2]))

  # Downsample
  if (sum(tumour$used) > N)
  {
    message(
      "\nMore than n = ",
      N,
      ' tumour mutations in the required VAF range, downsampling the data.\n'
    )

    w_t =  tumour %>% dplyr::filter(used) %>% dplyr::sample_n(N) %>% dplyr::pull(id)
    tumour$used[!(tumour$id %in% w_t)] = FALSE
  }

  t_ids = tumour %>% dplyr::filter(used) %>% dplyr::pull(id)

  germline$used = FALSE
  germline$used[germline$id %in% t_ids] = TRUE

  joint_data$used = FALSE
  joint_data$used[joint_data$id %in% t_ids] = TRUE

  # Print
  TB = table(tumour$used)
  pio::pioStr(
    "\n     VAF range",
    VAF_range_tumour[1],
    '~',
    VAF_range_tumour[2],
    ' [',
    paste(names(TB), TB),
    '] tumour mutations',
    suffix = '\n'
  )


  return(
    list(
      tumour = tumour,
      normal = germline,
      joint = joint_data,
      VAF_range_tumour = VAF_range_tumour,
      file = file
    )
  )
}

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
#' @export
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
    mutate(id = paste(chr, from, to, ref, alt, sep = ':'))

  # cli::cli_alert_info("Using for mutation data {.url {required_colnames}} for n = {.value {nrow(x)}}.")
  cli::cli_alert_success("Found data for n = {.value {nrow(x)}} mutations.")

  # Map CNA data, and retainn only mappable mutations
  if(!all(is.null(cna)))
  {
    cli::cli_alert_info("Found CNA data, mapping mutations to segmennts.")

    cn_obj = CNAqc::init(snvs = TINC:::as_tumour(x), cna = cna, .8)

    mappable = cn_obj$snvs %>%
      dplyr::filter(!is.na(segment_id)) %>%
      dplyr::pull(id)

    x = x %>%
      dplyr::filter(id %in% mappable)

    cli::cli_alert_success("Found {.value {nrow(cna)}} CNA segmennts, mapped n = {.value {nrow(x)}} mutations.")
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
    "Mutation with VAF within {.value {VAF_range_tumour}} ~ n = {.value { sum(x$OK_tumour) }}"
  )

  # Downsample
  if (sum(x$OK_tumour) > N)
  {
    cli::cli_alert_warning("More than {.value {N}} mutations, downsampling.")

    w_t =  x %>% dplyr::filter(OK_tumour) %>% dplyr::sample_n(N) %>% dplyr::pull(id)
    x$OK_tumour[!(x$id %in% w_t)] = FALSE
  }

  return(x)
}

as_normal = function(x) {
  x %>%
    dplyr::mutate(DP = n_ref_count + n_alt_count,
                  NV = n_alt_count,
                  VAF = NV / DP) %>%
    dplyr::select(-dplyr::ends_with('count'))
}

as_tumour = function(x) {
  x %>%
    dplyr::mutate(DP = t_ref_count + t_alt_count,
                  NV = t_alt_count,
                  VAF = NV / DP) %>%
    dplyr::select(-dplyr::ends_with('count'))
}

as_joint = function(x) {
  x %>%
    mutate(
      n_VAF = n_alt_count / (n_alt_count + n_ref_count),
      t_VAF = t_alt_count / (t_alt_count + t_ref_count)
    )
}

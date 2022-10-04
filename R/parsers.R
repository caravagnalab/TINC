#' Merged Canvas and Manta VCF parsing function
#'
#' @description Parse a VCF file merged from Canvas and
#' Manta, which stores a list of segments with absolute
#' CNA events, and tumour purity.
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
#' @importFrom stringr str_detect
#' 
#' @examples
#' # not run
#' \dontrun{
#'  load_VCF_Canvas_Manta("Myfile.vcf")
#' }
load_VCF_Canvas_Manta = function(file) {
  if (!file.exists(file))
    stop("Input file", file, "not found!")

  cli::cli_h2("VCF Canvas parser for TINC")

  # Read VCF with vcfR
  segments = vcfR::read.vcfR(file, verbose = FALSE)

  # Locations (Canvas only, not Manta)
  locs = segments@fix[, "ID"] %>%
    dplyr::as_tibble()
  row_ids = which(stringr::str_detect(locs$value, "Canvas"))
  locs = locs %>%
    filter(stringr::str_detect(value, "Canvas")) %>%
    tidyr::separate(value, into = c("A", "SD", "chr", "range"), sep = ":") %>%
    tidyr::separate(range, into = c("from", "to"), sep = "-") %>%
    dplyr::mutate(
      from = as.numeric(from),
      to = as.numeric(to)
      ) %>%
    select(-A, -SD)

  # Genotypes
  gt = vcfR::extract_gt_tidy(segments) %>%
    dplyr::mutate(
      Major = gt_MCC,
      minor = gt_CN - Major
    ) %>%
    select(minor, Major) %>%
    slice(row_ids)

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
#' \dontrun{
#'  load_VCF_Canvas("Myfile.vcf")
#' }
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



# Allele depth pull function
#
# @description Pulls out normal and tumour allele depths from
# Strelka 2 vcf columns
#
#
# @param CHROM as CHROM in vcf
# @param POS as POS in vcf
# @param REF as REF in vcf
# @param ALT as ALT in vcf
# @param FILTER as FILTER in vcf
# @param FORMAT as FORMAT in vcf
# @param normal as normal in vcf
# @param tumour as tumour in vcf
#
# @importFrom tibble tibble
# 
# @return Allele depths for normal and tumour
#
#
# @examples
# # not run
pullAD = function(CHROM, POS, REF, ALT, FILTER, FORMAT, normal, tumour){

#Identify ref and alt alleles
  refAlle = paste(REF, "U", sep = "")
  altAlle = paste(ALT, "U", sep = "")
  FORMAT_split = unlist(strsplit(FORMAT, ":"))
  refIndex = match(refAlle, FORMAT_split)
  altIndex = match(altAlle, FORMAT_split)

#Pull ref and alt depths for normal and tumour
  normal = unlist(strsplit(normal, ":"))
  tumour = unlist(strsplit(tumour, ":"))
  varAD = tibble::tibble(paste(CHROM,
                   as.numeric(POS) - 1,
                   POS,
                   REF,
                   ALT, sep = ":"),
                   strsplit(normal[refIndex], ",")[[1]][1],
                   strsplit(normal[altIndex], ",")[[1]][1],
                   strsplit(tumour[refIndex], ",")[[1]][1],
                   strsplit(tumour[altIndex], ",")[[1]][1],
                   FILTER
                   ) %>%
  purrr::set_names(c("id", "n_ref_count", "n_alt_count", "t_ref_count", "t_alt_count", "FILTERS")) %>%
  filter((as.numeric(t_alt_count) / (as.numeric(t_ref_count) + as.numeric(t_alt_count))) > 0,
         (as.numeric(t_alt_count) / (as.numeric(t_ref_count) + as.numeric(t_alt_count))) < 1)
 

  return(varAD)
}


#' Strelka 2 VCF parsing function
#'
#' @description Parse a VCF file from Srelka 2, which stores for
#' PASS SNVs the read depth of the alleles in the
#' nornmal and tumour
#'
#' @param file is Strelka 2 vcf file path
#'
#' @return Allele depths for normal and tumour
#' for all PASS autosome SNVs
#'
#' @export
#'
#' @importFrom  purrr  pmap_dfr
#'
#' @examples
#' # not run
#' \dontrun{
#'  load_VCF_Strelka("Myfile.vcf")
#' }
load_VCF_Strelka = function(file) {

#  file = "/home/jmitchell1/TIN/berthaRscript/input/LP3000417-DNA_F02_LP3000396-DNA_F02_12.vcf.gz"

  if (!file.exists(file))
    stop("Input file", file, "not found!")

# Load vcf and filter to leave PASS SNVs in autosome
  vcfSmallVar = read.table(file, colClasses = "character")
  autosome = sprintf("chr%s",seq(1:22))
  vcfSmallVarFilt = vcfSmallVar %>%
    dplyr::filter(V1 %in% autosome, nchar(V4) == 1, nchar(V5) == 1, V7 == "PASS")

# Pull allele depths from normal and tumour
  SNVnADtAD <- vcfSmallVarFilt[, c(1, 2, 4, 5, 7, 9, 10, 11)] %>%
    dplyr::rename_all(~ c("CHROM", "POS", "REF", "ALT", "FILTER", "FORMAT", "normal", "tumour")) %>%
    purrr::pmap_dfr(., pullAD) %>%
    separate(id, into = c('chr', 'from', 'to', 'ref', 'alt'), sep = ':') %>%
    mutate(from = as.numeric(from), to = as.numeric(to),
         n_ref_count = as.numeric(n_ref_count), n_alt_count = as.numeric(n_alt_count),
         t_ref_count = as.numeric(t_ref_count), t_alt_count = as.numeric(t_alt_count)) %>%
    select(-FILTERS)

  return(SNVnADtAD)
}


# Allele depth pull function for DRAGEN
#
# @description Pulls out normal and tumour allele depths from
# DRAGEN  vcf columns
#
#
# @param CHROM as CHROM in vcf
# @param POS as POS in vcf
# @param REF as REF in vcf
# @param ALT as ALT in vcf
# @param FILTER as FILTER in vcf
# @param FORMAT as FORMAT in vcf
# @param normal as normal in vcf
# @param tumour as tumour in vcf
#
# @importFrom tibble tibble
# 
# @return Allele depths for normal and tumour
#
#
# @examples
# # not run
pullAD_DRAGEN = function(CHROM, POS, REF, ALT, FILTER, FORMAT, normal, tumour){
  
  #Identify ref and alt alleles
  FORMAT_split = unlist(strsplit(FORMAT, ":"))
  ADIndex = match("AD", FORMAT_split)

  #Pull ref and alt depths for normal and tumour
  normal = unlist(strsplit(normal, ":"))
  tumour = unlist(strsplit(tumour, ":"))
  
  varAD = tibble::tibble(paste(CHROM,
                               as.numeric(POS) - 1,
                               POS,
                               REF,
                               ALT, sep = ":"),
                         strsplit(normal[ADIndex], ",")[[1]][1],
                         strsplit(normal[ADIndex], ",")[[1]][2],
                         strsplit(tumour[ADIndex], ",")[[1]][1],
                         strsplit(tumour[ADIndex], ",")[[1]][2],
                         FILTER
  ) %>%
    purrr::set_names(c("id", "n_ref_count", "n_alt_count", "t_ref_count", "t_alt_count", "FILTERS")) %>%
    filter((as.numeric(t_alt_count) / (as.numeric(t_ref_count) + as.numeric(t_alt_count))) > 0,
           (as.numeric(t_alt_count) / (as.numeric(t_ref_count) + as.numeric(t_alt_count))) < 1)
  
  
  return(varAD)
}


# Allele depth pull function for DRAGEN
#
# @description Pulls out CNVs values from
# DRAGEN  vcf columns
#
#
# @param CHROM as CHROM in vcf
# @param POS as POS in vcf
# @param REF as REF in vcf
# @param ALT as ALT in vcf
# @param FILTER as FILTER in vcf
# @param FORMAT as FORMAT in vcf
# @param sample as [normal_sample_name] in vcf
#
# @importFrom tibble tibble
# 
# @return CNA value parsed
#
#
# @examples
# # not run
pullAD_CNA_DRAGEN = function(chr, from, to, FORMAT, samples){
  
  #Identify ref and alt alleles
  FORMAT_split = unlist(strsplit(FORMAT, ":"))
  TotCNIndex = match("CN", FORMAT_split)
  MinCNIndex = match("MCN", FORMAT_split)
  
  
  #Pull ref and alt depths for normal and tumour
  samples = unlist(strsplit(samples, ":"))

  varAD = tibble::tibble(chr, from, to,
                         Major = as.numeric(strsplit(samples[TotCNIndex], ",")[[1]][1]),
                         minor= as.numeric(strsplit(samples[MinCNIndex], ",")[[1]][1])
  ) %>%  mutate(Major = Major - minor)
  
  
  return(varAD)
}



#' DRAGEN VCF small variants parsing function
#'
#' @description Parse a small variants VCF file from DRAGEN with information about SNVs and indels
#'
#' @param file is DRAGEN small variants vcf file path
#'
#' @return Allele depths for normal and tumour
#' for all PASS autosome SNVs
#'
#' @export
#'
#' @importFrom  purrr  pmap_dfr
#'
#' @examples
#' # not run
#' \dontrun{
#'  load_VCF_smallVar_DRAGEN("Myfile.vcf")
#' }
load_VCF_smallVar_DRAGEN = function(file) {
  

  if (!file.exists(file))
    stop("Input file", file, "not found!")
  
  cli::cli_h2("Small Variants VCF DRAGEN parser for TINC")
  
  # Load vcf and filter to leave PASS SNVs in autosome
  vcfSmallVar = read.table(file, colClasses = "character")
  autosome = sprintf("chr%s",seq(1:22))
  vcfSmallVarFilt = vcfSmallVar %>%
    dplyr::filter(V1 %in% autosome, nchar(V4) == 1, nchar(V5) == 1, V7 == "PASS")
  
  # Pull allele depths from normal and tumour
  SNVnADtAD <- vcfSmallVarFilt[, c(1, 2, 4, 5, 7, 9, 10, 11)] %>%
    dplyr::rename_all(~ c("CHROM", "POS", "REF", "ALT", "FILTER", "FORMAT", "normal", "tumour")) %>%
    purrr::pmap_dfr(., pullAD_DRAGEN) %>%
    separate(id, into = c('chr', 'from', 'to', 'ref', 'alt'), sep = ':') %>%
    mutate(from = as.numeric(from), to = as.numeric(to),
           n_ref_count = as.numeric(n_ref_count), n_alt_count = as.numeric(n_alt_count),
           t_ref_count = as.numeric(t_ref_count), t_alt_count = as.numeric(t_alt_count)) %>%
    select(-FILTERS)
  
  # Report some text
  cli::cli_alert_success(
    "DRAGEN variants parsed successfully: n = {.value {nrow(SNVnADtAD)}}"
  )
  
  
  return(SNVnADtAD)
}



#' DRAGEN VCF CNA parsing function
#'
#' @description Parse a VCF CNA file from DRAGEN, which has information on both 
#'
#' @param file is DRAGEN CNV vcf file path
#'
#' @return PASS copy number alteration values in CNAqc format, stored in the `cna`
#' slot, and estimated tumour purity in the `purity` slot.
#'
#' @export
#'
#' @importFrom  purrr  pmap_dfr
#'
#' @examples
#' # not run
#' \dontrun{
#'  load_VCF_CNA_DRAGEN("Myfile.vcf")
#' }
load_VCF_CNA_DRAGEN = function(file) {
  
  
  if (!file.exists(file))
    stop("Input file", file, "not found!")
  
  cli::cli_h2("CNV VCF DRAGEN parser for TINC")
  
  purity <- read.delim(file, sep = "\n") %>% grep(pattern = "##EstimatedTumorPurity", x = .[,1], value = T) %>%  strsplit(., split = "=")
  
  purity <- purity[[1]][2] %>%  as.numeric()
  
  # Load vcf and filter to leave PASS SNVs in autosome
  vcfSmallVar = read.table(file, colClasses = "character")
  vcfSmallVarFilt = vcfSmallVar %>%
    dplyr::filter(V7 == "PASS")  %>% 
    separate(V3, into = c("dragen", "type", "chr", "from_to"), sep = ":") %>% separate(from_to, into = c("from", "to"), sep = "-")
  
  # Pull allele depths from normal and tumour
  SNVnADtAD <- vcfSmallVarFilt %>% select(chr,from,to,V9,V10) %>%
    dplyr::rename_all(~ c("chr", "from", "to", "FORMAT", "sample")) %>%
    purrr::pmap_dfr(., pullAD_CNA_DRAGEN) 
  
  # Report some text
  cli::cli_alert_success(
    "DRAGEN parsed successfully: n = {.value {nrow(SNVnADtAD)}} CN segments."
  )
  
 
  cli::cli_alert_success(
    "DRAGEN tumour purity: p = {.value {purity}}."
  )
  
  
  return(list(cna = SNVnADtAD, purity = purity))
}


#' Check for CHIP variants in normal and tumour at high VAF
#'
#' @param normal_bam 
#' @param tumour_bam 
#' @param normal_vcf 
#' @param tumour_vcf 
#' @param reference 
#' @param plot 
#'
#' @return
#' @export
#'
#' @examples
check_CHIP_variants <- function(normal_bam = NULL, 
                                tumour_bam = NULL,
                                normal_vcf = NULL,
                                tumour_vcf = NULL,
                                reference = "hg38", 
                                plot = FALSE)
{
  
  if(all(is.null(normal_bam),is.null(tumour_bam), is.null(normal_vcf), is.null(tumour_vcf) )) stop("Please provide at least one of the inputs!")
  data('chip_hotspots')
  if(!(reference %in% c("hg38", "hg19")) stop("Please provide one of hg38 or hg19 as reference.")
  chip_variants <- chip_hotspots[[reference]]
  if(all(is.null(normal_bam), is.null(normal_vcf))){
    cli::cli_alert_danger("No normal vcf or bam provided, you will have no information abaout CHIP presence in the normal tissue!")
  } else if (!is.null(normal_bam)) {
    param <- Rsamtools::ScanBamParam(which=chip_variants)
    vars <- Rsamtools::pileup(normal_bam, scanBamParam = param)
  } else {
    
  }
  
}
#' Check for CHIP variants in normal and tumour at high VAF
#'
#' @param normal_bam path of bam file for healthy sample
#' @param tumour_bam path of bam file for tumour sample
#' @param reference reference genome, one of hg38 or hg19
#' @param min_vaf minimum required VAF to flag the mutation
#' @param only_snvs use only SVNs and filter indels 
#' @param add_chr if chromosomes on BAM files have chr set it to TRUE
#'
#' The function takes a tumour and a normal bam files and scans for known CHIP mutations listed in \href{https://www.nature.com/articles/s41586-020-2819-2}{Bick et al.}. The function gives back a data frame with mutation having VAF higher than \code{min_vaf}.
#' A bit of caution has to be taken when interpreting indels results, as from the pileup we cannot assign a lenght to the indel and the frameshit nature of the alteration. We suggest doing a proper small variant calling to confirm the results. 
#' @return Two dataframes for each of tumour and nomal with coloumns: chr, pos, alt, ref, VAF, DP and NV (chromosome, position, alternative allele, reference allele, variant allelic frequency, depth , number of variant reads)
#' @export
#'
check_CHIP_variants <- function(normal_bam = NULL, 
                                tumour_bam = NULL,
                                reference = "hg38", 
                                min_vaf = 0.05,
                                only_snvs = TRUE,
                                add_chr = TRUE)
{
  
  if(all(is.null(normal_bam),is.null(tumour_bam)))
    stop("Please provide at least one of the inputs!")
  data('chip_mutations')
  if(!(reference %in% c("hg38", "hg19"))) stop("Please provide one of hg38 or hg19 as reference.")
  chip_variants <- chip_mutations[[reference]]
  if(add_chr) chip_variants$chr = paste0("chr", chip_variants$chr)
  chip_variants_grange <- GenomicRanges::makeGRangesFromDataFrame(chip_variants, keep.extra.columns = T, seqnames.field = "chr", start.field = "pos", end.field = "pos")
  
  count_final_norm <- NULL
  if(is.null(normal_bam)){
    cli::cli_alert_danger("No normal bam provided, you will have no information abaout CHIP presence in the normal tissue!")
  } else {
    count_final_norm <-  check_CHIP_variants_BAM(normal_bam, chip_variants_grange, only_snvs, min_vaf)
  } 
  if(!is.null(count_final_norm)) cli::cli_alert_info("Found {nrow(count_final_norm)} CHIP associated variants in normal sample with VAF > {min_vaf}")
  
  count_final_tumour <- NULL
  if(is.null(tumour_bam)){
    cli::cli_alert_danger("No tumour bam provided, you will have no information abaout CHIP presence in the normal tissue!")
  } else {
    check_CHIP_variants_BAM(normal_bam, chip_variants_grange, only_snvs, min_vaf)
  } 
  if(!is.null(count_final_tumour)) cli::cli_alert_info("Found {nrow(count_final_tumour)} CHIP associated variants in tumour sample with VAF > {min_vaf}")
  
  
  return(list("normal" = count_final_norm, "tumour" = count_final_tumour))
}




check_CHIP_variants_BAM <- function(bam, chip_variants_grange, only_snvs, min_vaf){
  
  has_rsmatools <- require(Rsamtools)
  
  if(!has_rsmatools) stop("Please install Rsamtools to extract CHIP variants from bam files!")
  
  param <- Rsamtools::ScanBamParam(which=chip_variants_grange)
  pileup_params <- Rsamtools::PileupParam(distinguish_strands= FALSE, include_deletions= TRUE, include_insertions=
                                            TRUE)
  vars <- Rsamtools::pileup(bam, scanBamParam = param, pileupParam = pileup_params)
  count_table <- tidyr::pivot_wider(vars %>% dplyr::select(- which_label) %>% unique(), 
                                    id_cols = c("seqnames", "pos"), values_from = "count", 
                                    names_from = "nucleotide", values_fill = 0)
  count_table <- count_table %>% dplyr::rename(chr = seqnames,del = `-`, ins = `+` )
  count_merg <- dplyr::inner_join(count_table, chip_variants)
  count_merg <- count_merg %>% dplyr::mutate(type = dplyr::case_when(
    nchar(ref) > 1 ~ "del",
    nchar(alt) > 1 ~ "ins",
    TRUE ~ "SNP"
  ))
  count_merg_snps <- count_merg %>% filter(type == "SNP") 
  count_merg_snps$NV <-  sapply(1:nrow(count_merg_snps),function(x) count_merg_snps[x,count_merg_snps$alt[x]]) %>% unlist()
  count_merg_snps$DP <- count_merg_snps$NV + (sapply(1:nrow(count_merg_snps),function(x) count_merg_snps[x,count_merg_snps$ref[x]])  %>% unlist())
  count_merg_snps$VAF <- count_merg_snps$NV / count_merg_snps$DP
  count_final <- count_merg_snps %>% dplyr::filter(VAF > min_vaf)
  if(!only_snvs){
    count_merg_indel <- count_merg %>% filter(type != "SNP") 
    count_merg_indel$DP <-
      (sapply(1:nrow(count_merg_indel),function(x) sum(count_merg_indel[x,c("A", "C", "G", "T")]))  %>% unlist())
    
    count_merg_indel$NV <-  sapply(1:nrow(count_merg_indel),function(x) {
      if(count_merg_indel[x, "type"] == "del"){
        count_merg_indel[x,"del"]
      } else{
        count_merg_indel[x,"ins"]
      }
    }) %>% unlist()
    count_merg_indel$VAF <- count_merg_indel$NV / count_merg_indel$DP
    count_final_indel <- count_merg_indel %>% dplyr::filter(VAF > min_vaf)
    count_final <- rbind(count_final, count_final_indel)
  }
  
  return(count_final)
  
}



check_CHIP_variants_VCF <- function(vcf, chip_variants_grange, only_snvs, min_vaf){
  
  has_vcfR <- require(vcfR)
  
  if(!has_vcfR) stop("Please install vcfR to extract CHIP variants from vcf files!")
  
  vcf_read <- vcfR::read.vcfR(vcf)
  

  vcf_read_sub <- dplyr::inner_join(as.data.frame(vcf_read@fix) %>% mutate(POS = as.numeric(POS)),chip_variants %>% 
                                      dplyr::rename(CHROM = chr, POS = pos, ALT = alt, REF = ref)) 
  
  
    
}

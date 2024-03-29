% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsers.R
\name{load_VCF_smallVar_DRAGEN}
\alias{load_VCF_smallVar_DRAGEN}
\title{DRAGEN VCF small variants parsing function}
\usage{
load_VCF_smallVar_DRAGEN(file)
}
\arguments{
\item{file}{is DRAGEN small variants vcf file path}
}
\value{
Allele depths for normal and tumour
for all PASS autosome SNVs
}
\description{
Parse a small variants VCF file from DRAGEN with information about SNVs and indels
}
\examples{
# not run
\dontrun{
 load_VCF_smallVar_DRAGEN("Myfile.vcf")
}
}

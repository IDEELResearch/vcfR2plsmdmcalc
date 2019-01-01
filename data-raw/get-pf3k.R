#' P. falciparum 3,000 Genome Project
#'
#'
#' @format A vcf containing 4 samples and 100 biallelic SNPs
#' @section Warning: If you copy and paste the code below, it will write to
#'                   your Desktop/temp directory (if it exists).
#'
#' @section Dependencies: To generate this subsetted VCF, we used the
#'                        \code{NFBtools} package which is freely available
#'                        from GitHub with  \code{devtools::install_github("nickbrazeau/NFBtools")}
#'
#' @section Citation: This VCF has generously been made publicly available by
#' the MalariaGEN Consortium and the Pf3k project


url <- "ftp://ngs.sanger.ac.uk/production/pf3k/release_5/"
destfile <- "data-raw/="
httr::GET(url=url, httr::write_disk(path=destfile, overwrite = T))
pf3k <- vcfR::read.vcfR(file=destfile)
pf3k <- vcfR::extract.indels(pf3k[vcfR::is.biallelic(pf3k)], return.indels = F) # subset to SNPs










#' @title make_bifasta_mat
#'
#' @description makes a fasta matrix
#'
#' @param nsmpls number of strains to simulatie
#' @param nsites number of sites to simulate
#' @param npos_segsites A list of nucleotide positions to mutate
#' @param nseg_samples  A list of the number of samples to mutate at each of the positions. Note, lists are expected to be in same order
#'
#'


make_bifastamatvcfR <- function(nsmpls = 5, nsites = 100,
                             npos_segsites = list(),
                             nseg_samples = list()){
  #...........
  # error catch
  #...........
  if(length(unlist(npos_segsites)) <= 1){
  stop("You must provide at least two seg sites")
}
  if(max(unlist(npos_segsites)) > nsites){
    stop("The segregating site positision you supplied exceed haplotype length")
  }
  if(max(unlist(nseg_samples)) >= nsmpls){
    stop("The number of samples to mutate in your list exceeds the number minus 1, of samples you simulated")
  }
  #...........
  # simulate
  #...........
  dna <- c("a", "t", "c", "g")
  ref <- sample(dna, nsites, replace = T)
  dnamat <- matrix(rep(ref, each=nsmpls), nrow=nsmpls) # rows ind, col sites

  base <- sapply(unlist(npos_segsites), function(x){sample(
    dna[ which(! dna %in% dnamat[, x]) ], # pick another base
    size = 1)
    })

  for(i in 1:length(unlist(npos_segsites))){
    dnamat[ sample(1:nsmpls, size = nseg_samples[[i]], replace = F),
            npos_segsites[[i]]
            ] <- base[i] # anonymous for loop (not very good R)

  }

  fix <- data.frame(CHROM = "chrom",
             POS = unlist(npos_segsites),
             ID = NA,
             REF = toupper( ref[unlist(npos_segsites)] ),
             ALT = toupper( base ),
             QUAL = NA, FILT = NA, INFO = NA)

  compareDNA <- function(ref, muthap, npos_segsites){

    ret <- ifelse(ref[unlist(npos_segsites)] == muthap[unlist(npos_segsites)],
           "0", "1")
    return(ret)
  }

  gt <- apply(dnamat, 1, compareDNA, ref = ref, npos_segsites = npos_segsites)
  colnames(gt) <- paste0("smpl", seq(1:ncol(gt)))
  gt <- apply(gt, 2, function(x){paste0(x,":.:.")})
  gt <- cbind(data.frame(FORMAT = "GT:AD:DP"), gt)
  meta <- c("##fileformat=VCFv4.2", "##Simulated with vcfR2manip make_bifastamatvcfR")

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = as.matrix(fix), gt = as.matrix(gt))

  # make fasta from matrix
  dnamatlist <- apply(dnamat, 1, function(x) paste0(x, collapse = ""))

  names(dnamatlist) <- paste0("smpl", seq(1:length(dnamatlist)))
  dnamatlist <- Biostrings::DNAStringSet(dnamatlist)


  return(list(
    vcf = newvcfR,
    fasta = dnamatlist
  ))
}







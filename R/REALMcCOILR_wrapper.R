#' @title sample_posterior_mccoilr_categorical
#' @description Wrap the \code{REALMCCOILR} package functionality with posterior sampling

sample_posterior_mccoilr_categorical <- function(data, path = NULL, output = NULL,
                                                       totalrun=NULL, burnin=NULL,
                                                       ...){

  #................
  # error handle
  # ...............
  if(is.null(path)){
    stop("Must provide a path for the McCOIL_categorical function to write out to.")
  }
  if(is.null(output)){
    stop("Must provide a output for the McCOIL_categorical function to write out to.")
  }
  if(is.null(totalrun)){
    stop("Must provide a total run sampling for the McCOIL_categorical function.")
  }
  if(is.null(burnin)){
    stop("Must provide a burn-in for the McCOIL_categorical function.")
  }

  # run McCOIL_categorical
  z <- McCOILR::McCOIL_categorical(data = data, path=path, output=output, ... )
  # catch files to read in
  outputcatmoi <- list.files(path = path, pattern = output, full.names = T)
  outputcatmoi <- outputcatmoi[! grepl("_summary.txt", outputcatmoi)]
  if(purrr::is_empty(outputcatmoi)){
    stop(paste("The", output, "file was not created by
         the McCOILR::McCOIL_categorical function. Check your inputs"))
  }
  if(length(outputcatmoi) > 1){
    stop(paste("The", outputcatmoi, "file was found more than once. Duplicates not supported."))
  }

  #................
  # go
  # ...............
  nsampls <- nrow(data)
  outputcatmoi <- readr::read_tsv(outputcatmoi, col_names = F, comment = "total_acceptance")
  # file structure explanation here: https://github.com/Greenhouse-Lab/THEREALMcCOIL
  outputcatmoi <- outputcatmoi[(burnin+1):totalrun, 2:(nsampls+1)]
  colnames(outputcatmoi) <- rownames(data)
  ret <- apply(outputcatmoi, 2, function(x){ base::sample(x, size = 1, replace = F)})

  return(
    tibble::tibble(smpl = names(ret),
           MOI = ret)
  )

}


#' @title block_jackknife4sample_posterior_mccoilr_categorical
#' @description Wrap the \code{sample_posterior_mccoilr_categorical} functionality with posterior sampling for jackknife re-sampling

block_jackknife4sample_posterior_mccoilr_categorical <- function(data, block_size = 1, ...){

  ## nsites from McCOILR data is in columns
  nsites <- ncol(data)
  nsampls <- nrow(data)
  ## how many blocks, respecting (approx) the block size?
  M = (nsites/block_size)

  ## split into a list of arrays, each item a block
  lvls <- sort(rep(seq(1:floor(M)), block_size))
  lvls <- c(lvls, rep(NA, nsites-length(lvls))) # wasted ends
  blocks <- split(tibble::as.tibble(t(data)), lvls)
  M <- length(blocks)

  ## not enough blocks to get standard errors? too bad
  if(M < 3){

    warning("Not enough blocks to get standard errors. Returning NAs")

    return(list(
      se = NA,
      posterior_samples = NA
    ))
  }


  theta <- list()
  ## loop on blocks
  for(j in 1:length(blocks)){
    idx <- rep(T, length(blocks))
    idx[j] <- F
    keep <- t(dplyr::bind_rows(blocks[idx]))
    colnames(keep) <- paste0("V", 1:ncol(keep)) # mccoilr_categorical requires column names

    ## calculate statistic on remaining values
    theta[[j]] = sample_posterior_mccoilr_categorical(
        data = keep,
        output = paste0("rep", j, "-", "raw_McCOILR_output_catMOI.txt"),
        ...)
  }


  theta_df <- dplyr::bind_rows(theta, .id = "rep")

    ## calculate standard error
  theta_se <- theta_df %>%
    dplyr::group_by(smpl) %>%
    dplyr::mutate(theta_bar = mean(MOI)) %>%
    dplyr::summarise(sse = sum( (MOI - theta_bar)^2 ) )
  theta_se$se <- sqrt( ((M-1)/M) * theta_se$sse )

  ## return se and all posterior reps
  return(list(
    se = theta_se,
    posterior_samples = theta_df
  ))
}




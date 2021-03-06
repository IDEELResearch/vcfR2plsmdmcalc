#' @title block_jackknife
#' @description Block jackknife process from APM and vcfdo


block_jackknife <- function(values, statistic, block_size = 1, ...){

  if(! inherits(values, c("numeric", "integer")) ){
    stop(paste(class(values), " not supported for this operation.
               Must coerce to a (single) numeric or integer vector"))
  }

  ## how long is input? allow to be multidimensional, assume replicates are in rows
  n <- length(values)

  ## how many blocks, respecting (approx) the block size?
  M = (n/block_size)

  ## split into a list of arrays, each item a block
  lvls <- sort(rep(seq(1:floor(M)), block_size))
  lvls <- c(lvls, rep(NA, n-length(lvls))) # wasted ends
  blocks <- split(values, lvls)
  M <- length(blocks)

  ## not enough blocks to get standard errors? too bad
  if(M < 3){

    warning("Not enough blocks to get standard errors. Returning NAs")

    return(list(
      se = NA,
      theta = rep(NA, M)
    ))
  }

  theta <- rep(NA, M)


  ## loop on blocks
  for(j in 1:length(blocks)){
    idx <- rep(T, length(blocks))
    idx[j] <- F
    keep <- unlist(blocks[idx])
    ## calculate statistic on remaining values
    theta[j] <- statistic(keep)
  }

  ## calculate standard error
  theta_bar = mean(theta)
  se = sqrt( ((M-1)/M) * sum( (theta - theta_bar)^2 ) )
  #print("se = {}".format(se))

  ## return se and blockwise stats
  return(list(
    se = se,
    theta = theta
  ))
}

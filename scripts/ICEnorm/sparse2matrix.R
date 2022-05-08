sparse2matrix <- function(sparse, dim, resolution) 
{
  sparsetriangular <- sparse;
  if (ncol(sparsetriangular) >= 3) 
  {
    bins <- unique(c(sparsetriangular[, 1], sparsetriangular[, 2]))
    bins <- as.numeric(bins)
    bins <- bins[order(bins)]
    bins <- unique(bins)
    bin.size <- min(diff(bins))
    if(bin.size != resolution)
    {
      stop("Resolution assigned by user is not compatible with the corrodinates of bins")       
    }
    hicmatrix <- matrix(0, dim, dim)
    rownum <- sparsetriangular[, 1] / bin.size + 1
    colnum <- sparsetriangular[, 2] / bin.size + 1
    hicmatrix[cbind(rownum, colnum)] <- as.numeric(c(sparsetriangular[, 3]))
    hicmatrix[cbind(colnum, rownum)] <- as.numeric(c(sparsetriangular[, 3]))
    hicmatrix[is.na(hicmatrix)] <- 0
    
    return(hicmatrix)
  }
  else
  {
    stop("The number of columns of a sparse matrix should be not less than 3")       
  }
}
zerofilter <- function(mat) 
{
  zeros = unique(which(colSums(mat) == 0), which(rowSums(mat) == 0))
  if (length(zeros) > 0) 
  {
    rnames <- rownames(mat)[zeros]
    cnames <- colnames(mat)[zeros]
    mat = mat[-zeros, -zeros]
    filtedobj <- list(matrix=mat, rownames=rnames, colnames=cnames, lines=zeros)
  }
  else
    filtedobj <- list(matrix=mat, rownames=NULL, colnames=NULL, lines=NULL)
  return (filtedobj)
}
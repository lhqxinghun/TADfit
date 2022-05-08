funDownsample <- function(hicmat, n_sample){
  hicmat_rows <- nrow(hicmat)
  hicmat_cols <- ncol(hicmat)
  n = hicmat_rows*(hicmat_cols+1)/2
  
  hic_vec <- c()
  hicvec_pos <- c(1:n)
  #hicvec_x <- c()
  #hicvec_y <- c()
  
  for(i in 1:hicmat_cols)
  {
    hicmat_rows_tri <- i
    subvec <- hicmat[hicmat_rows_tri:hicmat_rows, i]
    hic_vec <- c(hic_vec, subvec)
  }
  
  hicvec_pos_rep = rep(hicvec_pos,hic_vec)
  lenrep = length(hicvec_pos_rep)
  
  x <- lenrep/n_sample
  
  
  sam_num <- sample(hicvec_pos_rep, x, replace=F)
  hicvec_sample <- 0*c(1:n)
  for (k in 1:lenrep)
  {
    sam_pos <- sam_num[k]
    hicvec_sample[sam_pos] <- hicvec_sample[sam_pos]+1
  }
  
  hicmat_sample <- matrix(0, hicmat_rows, hicmat_cols)
  for(i in 1:hicmat_cols)
  {
    hicmat_rows_tri <- i
    hicmat_rows_down <- (hicmat_rows+hicmat_rows-i+1)*i/2
    hicmat_rows_up <- hicmat_rows_down-hicmat_rows+i
    hicmat_sample[hicmat_rows_tri:hicmat_rows,i] <- hicvec_sample[hicmat_rows_up:hicmat_rows_down]
  }
  hicmat_sample[upper.tri(hicmat_sample)] <- t(hicmat_sample)[upper.tri(hicmat_sample)]
  return(hicmat_sample)
}


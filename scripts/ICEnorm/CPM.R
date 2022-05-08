CPM <- function(dirin, topfilenum, dirout, species, chr, resolution, bsparse, bcpm, bmodified) 
{
  switch(species,
    hg19=
    {
      chrlens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                  146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                  102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                  155270560, 59373566)
      names(chrlens) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                          "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                          "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
                          "chrX", "chrY")
      dim <- floor(chrlens[chr] / resolution) + 1      
    },
    mm9=
    {
      chrlens = c(197195432, 181748087, 159599783, 155630120, 152537259, 149517037,	152524553,
                  131738871,	124076172,	129993255,	121843856, 121257530,	120284312,	125194864,
                  103494974,	98319150,	95272651,	90772031,	61342430,	166650296,	15902555,	16300)
      names(chrlens) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                          "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                          "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY", "chrM")
      dim <- floor(chrlens[chr] / resolution) + 1      
    }    
  )
  datalist <- sapply(list.files(dirin, pattern=".txt", full.names=TRUE), read.table, sep = "\t", simplify = FALSE, USE.NAMES = TRUE)
  datalist <- datalist[1:topfilenum]
  if(bsparse == TRUE)
  {
    datalist <- lapply(datalist, sparse2matrix, species=species,chr=chr,resolution=resolution)
  }
  namesplit <- strsplit(basename(names(datalist)),split=".",fixed=TRUE)
  names(datalist) <- unlist(lapply(namesplit,head,1))
  CPMfun <- function(i)
  {
    mat <- as.matrix(datalist[[i]])
    mat[is.na(mat)] <- 0
    if(bcpm == TRUE)
    {
      if(bmodified == TRUE)
	{
	  libs <- sum(mat[upper.tri(mat, diag = TRUE)])
	  mat <- mat/libs*10^7
	}
	else
	{
	  libs <- sum(mat[upper.tri(mat, diag = TRUE)])
	  mat <- (mat+0.5)/(libs+1)*10^6
        }    
    }
    write.table(round(mat, digits=6), paste0(dirout, "/", basename(names(datalist)[i]), ".txt"), row.names = FALSE, col.names = FALSE, sep = "\t")
  }
  lapply(seq_along(datalist), CPMfun)
}

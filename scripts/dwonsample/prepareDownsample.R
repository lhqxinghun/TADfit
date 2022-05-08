rm(list = ls())
source("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/funDownsample.R")
source("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/submatrix.R")
source("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/sparse2matrix.R")

########################10K
#########GM12878
dirin <- "/media/disk2/dump/GM12878/chr1/10K"
n_downsample <- c(2, 4, 8, 16)
dim <- 24926 #for chr1 10K
filelist <- list.files(dirin, pattern = ".txt", full.names = TRUE)
for (file in filelist[1:5]) {
  hicmat <- as.matrix(sparse2matrix(read.table(file), dim, 10000))
  hicmat <- submatrix(hicmat, "hg19", "chr1", start = 150000000, end = 155000000, resolution = 10000)
  for (n in n_downsample) {
    dirout <- paste0("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/GM12878/chr1/10K-5M/10K-150M-155M/seq-depth_1_", n, "/Raw")
    hicmat_sample <- funDownsample(hicmat, n)
    fileout <- paste(dirout, basename(file), sep = "/")
    write.table(hicmat_sample, fileout, col.names = FALSE, row.names = FALSE, sep = "\t")
  }
}

#########IMR90
dirin <- "/media/disk2/dump/IMR90/chr1/10K"
n_downsample <- c(2, 4, 8, 16)
dim <- 24926 #for chr1 10K
filelist <- list.files(dirin, pattern = ".txt", full.names = TRUE)
for (file in filelist[1:5]) {
  hicmat <- as.matrix(sparse2matrix(read.table(file), dim, 10000))
  hicmat <- submatrix(hicmat, "hg19", "chr1", start = 150000000, end = 155000000, resolution = 10000)
  for (n in n_downsample) {
    dirout <- paste0("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/IMR90/chr1/10K-5M/10K-150M-155M/seq-depth_1_", n, "/Raw")
    hicmat_sample <- funDownsample(hicmat, n)
    fileout <- paste(dirout, basename(file), sep = "/")
    write.table(hicmat_sample, fileout, col.names = FALSE, row.names = FALSE, sep = "\t")
  }
}

#########K562
dirin <- "/media/disk2/dump/K562/chr1/10K"
n_downsample <- c(2, 4, 8, 16)
dim <- 24926 #for chr1 10K
filelist <- list.files(dirin, pattern = ".txt", full.names = TRUE)
for (file in filelist[1:5]) {
  hicmat <- as.matrix(sparse2matrix(read.table(file), dim, 10000))
  hicmat <- submatrix(hicmat, "hg19", "chr1", start = 150000000, end = 155000000, resolution = 10000)
  for (n in n_downsample) {
    dirout <- paste0("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/K562/chr1/10K-5M/10K-150M-155M/seq-depth_1_", n, "/Raw")
    hicmat_sample <- funDownsample(hicmat, n)
    fileout <- paste(dirout, basename(file), sep = "/")
    write.table(hicmat_sample, fileout, col.names = FALSE, row.names = FALSE, sep = "\t")
  }
}

rm(list = ls())
source("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/funDownsample.R")
source("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/submatrix.R")
source("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/sparse2matrix.R")
########################25K
#########GM12878
dirin <- "/media/disk2/dump/GM12878/chr1/25K"
n_downsample <- c(2, 4, 8, 16)
dim <- 9971 #for chr1 25K
filelist <- list.files(dirin, pattern = ".txt", full.names = TRUE)
for (file in filelist[1:5]) {
  hicmat <- as.matrix(sparse2matrix(read.table(file), dim, 25000))
  hicmat <- submatrix(hicmat, "hg19", "chr1", start = 150000000, end = 155000000, resolution = 25000)
  for (n in n_downsample) {
    dirout <- paste0("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/GM12878/chr1/25K-5M/25K-150M-155M/seq-depth_1_", n, "/Raw")
    hicmat_sample <- funDownsample(hicmat, n)
    fileout <- paste(dirout, basename(file), sep = "/")
    write.table(hicmat_sample, fileout, col.names = FALSE, row.names = FALSE, sep = "\t")
  }
}

#########IMR90
dirin <- "/media/disk2/dump/IMR90/chr1/25K"
n_downsample <- c(2, 4, 8, 16)
dim <- 9971 #for chr1 25K
filelist <- list.files(dirin, pattern = ".txt", full.names = TRUE)
for (file in filelist[1:5]) {
  hicmat <- as.matrix(sparse2matrix(read.table(file), dim, 25000))
  hicmat <- submatrix(hicmat, "hg19", "chr1", start = 150000000, end = 155000000, resolution = 25000)
  for (n in n_downsample) {
    dirout <- paste0("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/IMR90/chr1/25K-5M/25K-150M-155M/seq-depth_1_", n, "/Raw")
    hicmat_sample <- funDownsample(hicmat, n)
    fileout <- paste(dirout, basename(file), sep = "/")
    write.table(hicmat_sample, fileout, col.names = FALSE, row.names = FALSE, sep = "\t")
  }
}

#########K562
dirin <- "/media/disk2/dump/K562/chr1/25K"
n_downsample <- c(2, 4, 8, 16)
dim <- 9971 #for chr1 25K
filelist <- list.files(dirin, pattern = ".txt", full.names = TRUE)
for (file in filelist[1:5]) {
  hicmat <- as.matrix(sparse2matrix(read.table(file), dim, 25000))
  hicmat <- submatrix(hicmat, "hg19", "chr1", start = 150000000, end = 155000000, resolution = 25000)
  for (n in n_downsample) {
    dirout <- paste0("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/K562/chr1/25K-5M/25K-150M-155M/seq-depth_1_", n, "/Raw")
    hicmat_sample <- funDownsample(hicmat, n)
    fileout <- paste(dirout, basename(file), sep = "/")
    write.table(hicmat_sample, fileout, col.names = FALSE, row.names = FALSE, sep = "\t")
  }
}

rm(list = ls())
source("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/funDownsample.R")
source("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/submatrix.R")
source("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/sparse2matrix.R")
########################50K
#########GM12878
dirin <- "/media/disk2/dump/GM12878/chr1/50K"
n_downsample <- c(2, 4, 8, 16)
dim <- 4986 #for chr1 50K
filelist <- list.files(dirin, pattern = ".txt", full.names = TRUE)
for (file in filelist[1:5]) {
  hicmat <- as.matrix(sparse2matrix(read.table(file), dim, 50000))
  hicmat <- submatrix(hicmat, "hg19", "chr1", start = 150000000, end = 155000000, resolution = 50000)
  for (n in n_downsample) {
    dirout <- paste0("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/GM12878/chr1/50K-5M/50K-150M-155M/seq-depth_1_", n, "/Raw")
    hicmat_sample <- funDownsample(hicmat, n)
    fileout <- paste(dirout, basename(file), sep = "/")
    write.table(hicmat_sample, fileout, col.names = FALSE, row.names = FALSE, sep = "\t")
  }
}

#########IMR90
dirin <- "/media/disk2/dump/IMR90/chr1/50K"
n_downsample <- c(2, 4, 8, 16)
dim <- 4986 #for chr1 50K
filelist <- list.files(dirin, pattern = ".txt", full.names = TRUE)
for (file in filelist[1:5]) {
  hicmat <- as.matrix(sparse2matrix(read.table(file), dim, 50000))
  hicmat <- submatrix(hicmat, "hg19", "chr1", start = 150000000, end = 155000000, resolution = 50000)
  for (n in n_downsample) {
    dirout <- paste0("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/IMR90/chr1/50K-5M/50K-150M-155M/seq-depth_1_", n, "/Raw")
    hicmat_sample <- funDownsample(hicmat, n)
    fileout <- paste(dirout, basename(file), sep = "/")
    write.table(hicmat_sample, fileout, col.names = FALSE, row.names = FALSE, sep = "\t")
  }
}

#########K562
dirin <- "/media/disk2/dump/K562/chr1/50K"
n_downsample <- c(2, 4, 8, 16)
dim <- 4986 #for chr1 50K
filelist <- list.files(dirin, pattern = ".txt", full.names = TRUE)
for (file in filelist[1:5]) {
  hicmat <- as.matrix(sparse2matrix(read.table(file), dim, 50000))
  hicmat <- submatrix(hicmat, "hg19", "chr1", start = 150000000, end = 155000000, resolution = 50000)
  for (n in n_downsample) {
    dirout <- paste0("/media/disk1/liuerhu/hicda_liuerhu/seq-depth/K562/chr1/50K-5M/50K-150M-155M/seq-depth_1_", n, "/Raw")
    hicmat_sample <- funDownsample(hicmat, n)
    fileout <- paste(dirout, basename(file), sep = "/")
    write.table(hicmat_sample, fileout, col.names = FALSE, row.names = FALSE, sep = "\t")
  }
}

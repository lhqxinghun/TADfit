# TADfit

TADfit is a multivariate linear regression model for profiling of hierarchical chromatin domains, which tries to fit the interaction frequencies in Hi-C contact matrix with and without replicates using all-possible hierarchical TADs, and the regression coefficient for each TAD is obtained by an online learning solver called Follow-The-Regularized -Leader (FTRL) to determine the significant ones from them.

## Installation
 The TADfit was built on R 3.6.2 with gcc 7.4.0 and R 3.6.3 with gcc 8.1.0 for linux and windows respectively, and the same version of platforms is recommended for users if the compiled package is used.
 1) If necessary, install the required dependencies in R. 
```R
install.packages(c('Matrix', 
                'perm', 
                'matrixStats', 
                'Rcpp', 
                'HiCcompare', 
                'plyr', 
                'zoo'))
```
 2) Install TADfit from source package
```bash
  install.packages("path to TADfit_1.2.tar.gz", repos = NULL, type = "source") 
```  
  
 3) Or install TADfit from compiled package
```bash
#linux
install.packages('path to TADfit_1.2_R_x86_64-pc-linux-gnu.tar.gz')
```
```bash
#Windows
install.packages('path to TADfit_1.2.zip', repos = NULL, type = "win.binary")
```

## Usage
The hierarchical TADs on replicate Hi-C data are called mainly by the  function *calltad* in TADfit. Other auxiliary functions can be accessable in the package too.
### Input data
Before runing the TADfit, Users need to prepare the input data. The input data to TADfit is a group of replicate normalized Hi-C contact matrices under the same scale. Here we recommend iterative correction and eigenvector decomposition (ICE) followed by Counts per million (CPM) for normalization. These matrix files need to be put on the same directory.

### Parameters
- **dirin**: path to a directory which contains replicate Hi-C contact matrix file(s), Logarithmic operations is integrated in the althrithm and the input should be raw contact matrices or those normalized with ICE or other normalization methods followed by CPM (recommended)
- **resolution**: bin size of the Hi-C contact matrix
- **start**: start coordinate of the chromosome region
- **end**: end coordinate of the chromosome region
- **specied**: "hg19" for human or "mm9"for mouse
- **chr**: chromosome number, (default "chr1"), chr and species are used to checke the start and end coordinates
- **winsize**:window size when calling TAD boundaries. This parameter is from TopDom (default 5)
- **minsize**: the minimum size of allowed TADs measured in bins(default 3)
- **maxsize**: the maximum size of allowed TADs measured in bins(default 200)
- **repeatnum**: the number of repeat when solving the model(default 5)
- **iteration**: the number of iteration when solving the model (default 5 and can be tunned down to reduce time complexity but not recommended to be smaller than 2)
- **alpha, beta**: parameters for pre-coordinating learning rate when optimizing the model (default 0.01 for alpha and 1 for beta)
- **l**: adjust the strength of l1 regularization (defalut 2)

### Output
The model(*calltad* function) outputs a dataframe with four colomns
- **start**: start coordinate (in bin) of a hierarchical TAD
- **end**: end coordinate (in bin)of a hierarchical TAD
- **b_avg**: the average regression coefficients for each hierachical TAD estimated by FTRL
- **pvalue**：the p-value for permutation test

An example of output:
```
start   end     b_avg     pvalue
0       7       0.11833   0.00082
2       7       0.14427   0.00079
7       18      0.22934   0.00054
7       33      0.03613   0.00157
...
```

### Example
Call hierarchical TADs from the sample Hi-C data
```R
library(TADfit)
dirin <- "/media/disk1/liuerhu/hicda_liuerhu/github/data/expdata/IMR90-chr21-25K"
resolution <- 25000
start <- 0
end <- 48000000
TADs <- calltad(dirin, resolution, start, end)  #Other parameters are default
```
## Reference
Liu, E., Lyu, H., Peng, Q. et al. TADfit is a multivariate linear regression model for profiling hierarchical chromatin domains on replicate Hi-C data. Commun Biol 5, 608 (2022).

## Contact
hongqianglv@xjtu.edu.cn
liuerhu@stu.xjtu.edu.cn

## Update
Compared to TADfit1.1, the efficiency of this TADfit1.2 has been greatly improved by defaultly limiting the maximum size of the TADs to 200 bins, minimum size of the TADs to 3 bins, and removing redundant calculations without reducing the accuracy of the model.

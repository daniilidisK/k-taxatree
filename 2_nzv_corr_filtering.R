# Clear
cat("\014")
rm(list = ls())

# libraries
library(caret)
library(stats)

# fucntion - PCA feature selection -  NOT USED
PCA_feature_selection <- function(kmerMatrix, info_perc){
  
  kmerMatrix.pca <- prcomp(kmerMatrix, center = FALSE, scale. = FALSE)
  rot <- as.data.frame(kmerMatrix.pca$rotation)
  
  # keep features with cumulative PVE <= information percentage
  std_dev <- kmerMatrix.pca$sdev
  pr_var <- std_dev^2
  
  prop_varex <- pr_var/sum(pr_var)
  n_kmers <- 1
  count_perc <- 0
  
  while(count_perc <= info_perc) {
    count_perc <- count_perc + prop_varex[n_kmers]
    n_kmers <- n_kmers +1 
  }
  
  rot <- rot[, 1:(n_kmers-1)]
  
  for (i in 1:ncol(rot)){
    pc <- abs(rot[,i])
    mean_pc <- mean(pc)
    var_pc <- sqrt(var(pc))
    
    indices <- which(pc < mean_pc + 3*var_pc)
    rot[indices, i] <- 0
    
  }
  
  kmers_to_keep <- names(which(rowSums(rot) != 0))
  columns_to_keep <- which(colnames(kmerMatrix) %in% kmers_to_keep)
  kmerMatrix <- kmerMatrix[, columns_to_keep]
  
  return(kmerMatrix)
}

# Importing kmerMatrix
#kmerMatrix <- read.csv('C:/Users/User/Desktop/INAB/k-mers/empdata/emp-train-test_matrix_k_7.csv', 
#                       row.names = 1)

kmerMatrix <- read.csv('emp-data-loc/kmerMatrix.csv', row.names = 1)

# remove near zero variance predictors
nzv <- nearZeroVar(kmerMatrix, freqCut = 15, uniqueCut = 0.5) # 0.3 

if(length(nzv) > 0){
  kmerMatrix <- kmerMatrix[, -nzv]
  rm(nzv)
}

# Correlation
corMat <- cor(kmerMatrix)
highlyCorNZV <- findCorrelation(corMat, cutoff = .8) #0.85

# keep low correlated features
kmerMatrix <- kmerMatrix[, -highlyCorNZV]
rm(corMat, highlyCorNZV)

# PCA feature selection
kmerMatrix <- PCA_feature_selection(kmerMatrix, info_perc = 0.90)

dir.create("emp-data-loc")
# save
write.csv(kmerMatrix, 'emp-data-loc/kmerMatrix_nzv_corr_PCA_099.csv')

# Clear
cat("\014")
rm(list = ls())

# libraries
library(kmer)
library(seqinr)
library(stringr)
library(parallel)
library(gtools)
library(dplyr)
library(data.table)
library(stats)
library(cluster)
library(parallelDist)
library(RcppParallel)
library(microbenchmark)
library(doParallel)

detectedCores <- parallel::detectCores()
registerDoParallel(cores=detectedCores-1) 

######################### FUNTIONS #############################################
count_kmers_in_seq <- function(seq, k, seq_name){
  
  one.run <- function(i) {
    substr(seq, start = i, stop = i+k - 1)
  }
  
  out <- lapply(c(1:(nchar(seq)-k+1)), one.run)
  kmers_row <- plyr::count(unlist(out))
  kmers_row_mat <- matrix(kmers_row$freq, 
                          nrow = 1, 
                          ncol = nrow(kmers_row), 
                          dimnames = list(c(seq_name), kmers_row$x))
  kmers_row_mat <- as.data.frame(kmers_row_mat)
  return(kmers_row_mat)
}


count_kmers_in_file <- function(sequences, k, seq_names){
  
  count.kmers <- function(i){
    count_kmers_in_seq(sequences[i], k, seq_names[i])
  }
  out <- mclapply(c(1:length(sequences)), count.kmers, mc.cores = 4)
  
  kmerMatrix <- as.matrix(rbindlist(out, fill = T))
  kmerMatrix[which(is.na(kmerMatrix), arr.ind = T)] <- 0
  rownames(kmerMatrix) <- seq_names
  
  return(kmerMatrix)
}

pairwise_distances <- function(X){
  
  minkowski_dist <- function(a,b){
    
    one.run <- function(i){
      abs(a[i] - b[i])
    }
    
    parallelSum <- function(x) {
      
      items <- length(x)
      batches <- detectedCores * 4
      batchSets <- split(x, rep(1:batches, length.out=items))
      
      finalSum <- foreach(b=iter(batchSets, by='row'), .combine="+") %dopar% sum(b)
      
      return (finalSum)
    }
    
    #return(parallelSum(unlist(mclapply(c(1:length(a)), one.run, mc.cores = 4))))
    return(parallelSum(abs(a-b)))
    
    
    #sum(unlist(mclapply(c(1:length(a)), one.run, mc.cores = 4)))
    
  }
  
  calculate_distances <- function(X, pos){
    
    n <- nrow(X)
    this_row <- X[pos, ]
    
    one.run <- function(i){
      minkowski_dist(this_row, X[i, ])
    }
    
    return(unlist(lapply(c((pos+1): n-1), one.run )))
    
  }
  
  
  D <- matrix(0L, nrow = nrow(X), ncol = nrow(X))
  dim_D <- nrow(X)
  
  for (j in 1:(dim_D - 1)){
    print(paste0('Row ', j))
    row_dists <- calculate_distances(X, j)
    D[(j+1):dim_D, j] <- row_dists
    D[j,(j+1):dim_D] <- row_dists
  }
  
  return(D)
  
}


###################### ANALYSIS ################################################

# input files
taxa_file <- read.csv('emp-data-loc/emp-taxonomy-train-test-loc.csv')

# taxa
taxa <- c('kingdom', 'phylum', 'class', 'order')

# values of k
kvals <- c(4:8)

# sequences, seqnames
sequences <- taxa_file$sequence
seq_names <- taxa_file$ID

# initialization
silh_scores_matrix <- matrix(0L, nrow = length(kvals), ncol = length(taxa))
rownames(silh_scores_matrix) <- kvals
colnames(silh_scores_matrix) <- taxa

# Analysis
for (k in kvals){
  
  print(paste0("Analysis for k = ", k))
  kmerMatrix <- count_kmers_in_file(sequences = sequences, seq_names = seq_names, k = k)
  print('Kmer matrix created')
  
  D <- as.matrix(dist(kmerMatrix, method = 'minkowski', p = 1))
  
  print('Distance matrix calculated')
  rm(kmerMatrix)
  
  Sys.sleep(5)
  
  print('go')
  silh_row <- c()
  
  for (tax in taxa){
    
    print(paste0('Taxa: ', tax))
    
    targets <- taxa_file[,tax]
    to_drop <- which(targets == 'Unassigned')
    
    if(length(to_drop)>0){
      
      t <- factor(targets[-to_drop])
      s <- silhouette(as.numeric(t), D[-to_drop, -to_drop])
      silh_row <- c(silh_row, mean(s[,3]))
      
    } else {
      
      t <- factor(targets)
      s <- silhouette(as.numeric(t), D)
      silh_row <- c(silh_row, mean(s[,3]))
    
    }
    
  }
  
  silh_scores_matrix[as.character(k),] <- silh_row
  
  rm(silh_row, targets, D)
  
}

write.csv(silh_scores_matrix, file = 'silh.csv')


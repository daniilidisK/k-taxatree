# Clear
cat("\014")
rm(list = ls())

# libraries
library(stringr)
library(parallel)
library(gtools)
library(dplyr)
library(data.table)
library(stats)
library(cluster)
# library(reticulate)
library(stringdist)
library(Rfast)

# use_python("/home/togkousa/Anaconda3/envs/r-reticulate/python.exe")

# for server
# use_python("/home/user/anaconda3/envs/togkou_1/bin/python.exe")

# source_python('1_Kselection-tool/dist.py')


# detectedCores <- parallel::detectCores()
# registerDoParallel(cores=detectedCores-1) 

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
  out <- mclapply(c(1:length(sequences)), count.kmers, mc.cores = 8)
  
  kmerMatrix <- as.matrix(rbindlist(out, fill = T))
  kmerMatrix[which(is.na(kmerMatrix), arr.ind = T)] <- 0
  rownames(kmerMatrix) <- seq_names
  
  return(kmerMatrix)
}

pairwise_distances <- function(X){
  
  convert_to_char <- function(X){
    
    one.run <- function(i){
      paste(as.character(X[i,]), collapse = '')
    }
    
    unlist(mclapply(c(1:nrow(X)), one.run, mc.cores = 4 ))
    
  }
  
  s <- convert_to_char(X)
  
  D <- stringdistmatrix(s, s, 'hamming', nthread = 4)
  rownames(D) <- rownames(X)
  colnames(D) <- rownames(X)
  
  return(D)
}


###################### ANALYSIS ################################################

# input files
taxa_file <- read.csv('emp-data-loc/emp-taxonomy-train-test-loc.csv')

# taxa
taxa <- c('kingdom', 'phylum', 'class', 'order')

# values of k
kvals <- c(4:10)

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
  
  # calculate pairwise distances
  # D <-  as.matrix(dist(kmerMatrix, method = 'minkowski', p = 1))
  # D <- pairwise_distances(kmerMatrix)
  D <- Dist(kmerMatrix, method = 'minkowski', p = 1)
  
  print('Distance matrix calculated')
  rm(kmerMatrix)
  
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


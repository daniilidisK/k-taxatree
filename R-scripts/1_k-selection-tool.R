# Clear
cat("\014")
rm(list = ls())

# libraries
library(parallel)
library(data.table)
library(cluster)
library(Rfast)
library(plyr)

# INPUTS -----------------

# input files
taxa_file <- read.csv('emp-data/emp-taxonomy-train-test.csv')

# taxa
taxa <- c('kingdom', 'phylum', 'class', 'order')

# values of k
kvals <- c(9:12)

# output folder
dir.create("Output")
output_folder <- "Output/k-selection"
dir.create(output_folder)

# FUNTIONS --------------------

# count_kmers_in_seq <- function(seq, k, seq_name){
#   
#   one.run <- function(i) {
#     substr(seq, start = i, stop = i+k - 1)
#   }
#   
#   out <- lapply(c(1:(nchar(seq)-k+1)), one.run)
#   kmers_row <- count(unlist(out))
#   kmers_row_mat <- matrix(kmers_row$freq, 
#                           nrow = 1, 
#                           ncol = nrow(kmers_row), 
#                           dimnames = list(c(seq_name), kmers_row$x))
#   kmers_row_mat <- as.data.frame(kmers_row_mat)
#   return(kmers_row_mat)
# }


count_kmers_in_file <- function(sequences, k, seq_names){
  
  # count.kmers <- function(i){
  #   count_kmers_in_seq(sequences[i], k, seq_names[i])
  # }
  # out <- mclapply(c(1:length(sequences)), count.kmers, mc.cores = 8)
  # out <- lapply(c(1:length(sequences)), count.kmers)
  
  one.run <- function(i, seq, k) {
    substr(seq, start = i, stop = i+k - 1)
  }
  
  out <- list()
  
  for (one.seq in 1:length(sequences)) {
    
    temp <- lapply(c(1:(nchar(sequences[1])-k+1)), 
                   one.run, 
                   seq = sequences[one.seq], k = k)
    kmers_row <- count(unlist(temp))
    kmers_row <- matrix(kmers_row$freq, 
                        nrow = 1, 
                        ncol = nrow(kmers_row), 
                        dimnames = list(c(seq_names[one.seq]), kmers_row$x))
    
    out[[one.seq]] <- as.data.table(kmers_row)
  }
 
  one_fifth <- round(length(sequences) / 10)

  kmerMatrix <- data.table()

  for (one.part in 1:10){
    
    temp <- rbindlist(out[1:one_fifth], fill = T)
    
    # print(length(which(is.na(temp))))
    temp[is.na(temp)] <- 0
    out <- out[-c(1:one_fifth)]
    kmerMatrix <- rbind(kmerMatrix, temp, fill = T)
    
    # print(length(which(is.na(kmerMatrix))))
    kmerMatrix[is.na(kmerMatrix)] <- 0
    rm(temp)
    print(one.part)
  }

  if (length(out) > 0) {
    temp <- rbindlist(out[one_fifth *10: length(out)], fill = T)
    temp[is.na(temp)] <- 0
    kmerMatrix <- rbind(kmerMatrix, temp, fill = T)
  }
  
  rm(out, temp)
  
  # print(length(which(is.na(kmerMatrix))))

  # kmerMatrix <- as.matrix(kmerMatrix)
  
  # kmerMatrix <- as.matrix(rbindlist(out, fill = T))

  # kmerMatrix[which(is.na(kmerMatrix), arr.ind = T)] <- 0
  
  rownames(kmerMatrix) <- seq_names
  
  return(kmerMatrix)
}


#ANALYSIS ------------------------

# sequences, seqnames
sequences <- taxa_file$sequence
seq_names <- taxa_file$ID

# initialization
silh_scores_matrix <- matrix(0L, nrow = length(kvals), ncol = length(taxa))
rownames(silh_scores_matrix) <- kvals
colnames(silh_scores_matrix) <- taxa

output_txt <- paste0(output_folder, "/silhouette_scores.txt")
file.create(output_txt)
write(paste(taxa, collapse =","), output_txt, append = T)

# Analysis
for (k in kvals){
  
  print(paste0("Analysis for k = ", k))
  kmerMatrix <- count_kmers_in_file(sequences = sequences, seq_names = seq_names, k = k)
  print('Kmer matrix created')
  
  # calculate pairwise distances
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
  write(silh_row, output_txt, append = T)
  
  rm(silh_row, targets, D, t, s, to_drop)
  
}



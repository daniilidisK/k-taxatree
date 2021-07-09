
count_kmers_in_seq <- function(one.seq, sequences, k, seq_names ){
  
  one.run <- function(i, seq, k) {
    substr(seq, start = i, stop = i+k - 1)
  }
  
  temp <- lapply(c(1:(nchar(sequences[1])-k+1)), 
                 one.run, 
                 seq = sequences[one.seq], k = k)
  kmers_row <- count(unlist(temp))
  rm(temp)

  kmers_row <- matrix(kmers_row$freq,
                      nrow = 1,
                      ncol = nrow(kmers_row),
                      dimnames = list(c(seq_names[one.seq]), kmers_row$x))

  kmers_row <- as.data.table(kmers_row)
  
  return(kmers_row)
}


count_kmers_in_file <- function(sequences, k, seq_names){
  
  # out <- mclapply(c(1:length(sequences)),
  #                 count_kmers_in_seq, 
  #                 mc.cores = 8, 
  #                 sequences = sequences, k = k, seq_names = seq_names)
  
  out <- lapply(c(1:length(sequences)),
                count_kmers_in_seq,
                sequences = sequences, k = k, seq_names = seq_names)
  
  one_fifth <- round(length(sequences) / 20)
  
  kmerMatrix <- data.table()
  
  for (one.part in 1:20){
    
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
    temp <- rbindlist(out[one_fifth * 20: length(out)], fill = T)
    temp[is.na(temp)] <- 0
    kmerMatrix <- rbind(kmerMatrix, temp, fill = T)
    rm(temp)
  }
  
  rm(out)
  
  # print(length(which(is.na(kmerMatrix))))
  
  # kmerMatrix <- as.matrix(kmerMatrix)
  
  # kmerMatrix <- as.matrix(rbindlist(out, fill = T))
  
  # kmerMatrix[which(is.na(kmerMatrix), arr.ind = T)] <- 0
  
  rownames(kmerMatrix) <- seq_names
  
  return(kmerMatrix)
}

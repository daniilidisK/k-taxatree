# Clear
cat("\014")
rm(list = ls())

# libraries
library(stringdist)
library(Rfast)
library(reticulate)
library(parallel)
library(png)

use_python("/home/togkousa/Anaconda3/envs/r-reticulate/python.exe")
source_python('1_Kselection-tool/dist.py')

# string dist function
string_dist_R <- function(X){
  
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


############# Constant cols, Rows changing #####################################
dim_1 <- c(100, 200, 300, 500, 1000, 2000)
dim_2 <- 1000


distances_overall <- matrix(0L, nrow = length(dim_1), ncol = 4)
colnames(distances_overall) <- c("dist", "stringdist", "Dist-Rfast", "python_dist")
rownames(distances_overall) <- dim_1


for (this_dim in dim_1){
  
  print(paste0('Analysis for nrows =  ', this_dim))
  
  this_row <- c()
  
  # create data matrix
  data <- matrix(sample.int(3, size = this_dim*dim_2, replace = TRUE)-1, nrow = this_dim, ncol = dim_2)
  
  # dist
  start_time <- Sys.time()
  print('dist')
  d <- as.matrix(dist(data, method = 'minkowski', p = 1))
  end_time <- Sys.time()
  diff_time <- as.numeric(end_time - start_time)
  print(diff_time)
  this_row <- c(this_row, diff_time)
  rm(d)
  
  # stringdist
  print('stringdist')
  start_time <- Sys.time()
  d <- string_dist_R(data)
  end_time <- Sys.time()
  diff_time <- as.numeric(end_time - start_time)
  print(diff_time)
  this_row <- c(this_row, diff_time)
  rm(d)
  
  # Dist
  print('fast dist')
  start_time <- Sys.time()
  d <- Dist(data, method = 'minkowski', p = 1)
  end_time <- Sys.time()
  diff_time <- as.numeric(end_time - start_time)
  print(diff_time)
  this_row <- c(this_row, diff_time)
  rm(d)
  
  # python
  start_time <- Sys.time()
  print('python dist')
  d <- dist_python(data)
  end_time <- Sys.time()
  diff_time <- as.numeric(end_time - start_time)
  print(diff_time)
  this_row <- c(this_row, diff_time)
  rm(d)
  
  distances_overall[as.character(this_dim), ] <- this_row
}

# plot
png('rows_changing.png', width = 1200, height = 800)
plot(dim_1, distances_overall[,1], type="l",col="red", ylim = c(0, max(distances_overall) + 0.5), 
     xlab = 'Number of rows', ylab = 'execution time')
lines(dim_1, distances_overall[,2], type="l",col="green")
lines(dim_1, distances_overall[,3], type="l",col="blue")
lines(dim_1, distances_overall[,4], type="l",col="yellow")
title('Execution time comparison - constant num of columns, num of rows changing ')
legend('topright', legend = c('dist', 'stringdist', 'fast_dist', 'python_dist'), col = c('red','green', 'blue', 'yellow'), lty = 1)
dev.off()


############# Constant rows, Cols changing #####################################
dim_2 <- c(100, 200, 300, 500, 1000, 2000)
dim_1 <- 1000


distances_overall <- matrix(0L, nrow = length(dim_2), ncol = 4)
colnames(distances_overall) <- c("dist", "stringdist", "Dist-Rfast", "python_dist")
rownames(distances_overall) <- dim_2


for (this_dim in dim_2){
  
  print(paste0('Analysis for ncols =  ', this_dim))
  
  this_row <- c()
  
  # create data matrix
  data <- matrix(sample.int(3, size = dim_1*this_dim, replace = TRUE)-1, nrow = dim_1, ncol = this_dim)
  
  # dist
  start_time <- Sys.time()
  print('dist')
  d <- as.matrix(dist(data, method = 'minkowski', p = 1))
  end_time <- Sys.time()
  diff_time <- as.numeric(end_time - start_time)
  print(diff_time)
  this_row <- c(this_row, diff_time)
  rm(d)
  
  # stringdist
  print('stringdist')
  start_time <- Sys.time()
  d <- string_dist_R(data)
  end_time <- Sys.time()
  diff_time <- as.numeric(end_time - start_time)
  print(diff_time)
  this_row <- c(this_row, diff_time)
  rm(d)
  
  # Dist
  print('fast dist')
  start_time <- Sys.time()
  d <- Dist(data, method = 'minkowski', p = 1)
  end_time <- Sys.time()
  diff_time <- as.numeric(end_time - start_time)
  print(diff_time)
  this_row <- c(this_row, diff_time)
  rm(d)
  
  # python
  start_time <- Sys.time()
  print('python dist')
  d <- dist_python(data)
  end_time <- Sys.time()
  diff_time <- as.numeric(end_time - start_time)
  print(diff_time)
  this_row <- c(this_row, diff_time)
  rm(d)
  
  distances_overall[as.character(this_dim), ] <- this_row
}

# plot
png('cols_changing.png', width = 1200, height = 800)
plot(dim_2, distances_overall[,1], type="l",col="red", ylim = c(0, max(distances_overall) + 0.5), 
     xlab = 'Number of cols', ylab = 'execution time')
lines(dim_2, distances_overall[,2], type="l",col="green")
lines(dim_2, distances_overall[,3], type="l",col="blue")
lines(dim_2, distances_overall[,4], type="l",col="yellow")
title('Execution time comparison - constant num of rows, num of cols changing ')
legend('topright', legend = c('dist', 'stringdist', 'fast_dist', 'python_dist'), col = c('red','green', 'blue', 'yellow'), lty = 1)
dev.off()
# Clear
cat("\014")
rm(list = ls())

# libraries
library(parallel)
library(data.table)
library(plyr)

source("R-scripts/count_kmers_functions.R")

# INPUTS -----------------

# input files
taxa_file <- read.csv('emp-data/emp-taxonomy-train-test.csv')

# value of k
k <- 8

# output folder
dir.create("Output")
output_folder <- "Output/kmer-matrix-creation"
dir.create(output_folder)

# MATRIX CREATION ------------------------

kmerMatrix <- count_kmers_in_file(sequences = taxa_file$sequence, 
                                  seq_names = taxa_file$ID, 
                                  k = k)

fwrite(kmerMatrix, paste0(output_folder, "/kmerMatrix_k_8.csv"))
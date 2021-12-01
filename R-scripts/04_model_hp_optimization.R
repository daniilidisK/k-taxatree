# Clear
cat("\014")
rm(list = ls())
# dev.off(dev.list()["RStudioGD"])

# libraries
library(caret)
library(UBL)
library(splitstackshape)
library(mlr)
library(data.table)
library(mldr)
library(dplyr)
library(parallel)
library(hash)

# source
source('multilabel_functions.R')

# get input from user
# args = commandArgs(trailingOnly=TRUE)
mtry_list <- c(16, 17, 18, 19, 20)
ntrees_list <- seq(50, 500, 50)
num_of_experiments <- 10
predict_theshold <- 0.5
output_folder <- 'Output/hp optimization'

# taxa
taxa <- c('kingdom', 'phylum', 'class', 'order')

# put input data folder
input_folder <- "emp-data"

# initialization of output folders
previous_folder <- "Output/preprocessing"
dir.create(output_folder)

for (tax in taxa){
  to_create <- paste0(output_folder, "/", tax)
  dir.create(to_create)
}

# taxonomies table
taxonomyFilepath <- paste0(input_folder, '/emp-taxonomy-train-test.csv')
taxonomies_table <- read.csv(taxonomyFilepath) #[,-1]
# taxonomies_table <- taxonomies_table[1:400, ]

# Read data matrix
kmerMatrix <- read.csv(paste0(previous_folder, '/kmerMatrix_fs.csv'),
                       row.names = 1)
rownames(kmerMatrix) <- taxonomies_table$ID

# kmerMatrix <- kmerMatrix[1:400, ]

# Dropping totally Unassigned/ Unclassified data
out <- excluding_unassigned(kmerMatrix, taxonomies_table, taxa)
kmerMatrix <- out[[1]]
taxonomies_table <- out[[2]]

# keeping only kingdom, phylum, class, order
taxonomies_table <- taxonomies_table[,taxa]
rownames(taxonomies_table) <- rownames(kmerMatrix)

# Setting 'Unassigned' as ''
taxonomies_table[which(taxonomies_table == 'Unassigned', arr.ind = TRUE)] <- ''

# excluding singletons
out <- excluding_singletons(kmerMatrix, taxonomies_table, taxa)
kmerMatrix <- out[[1]]
taxonomies_table <- out[[2]]

# oversampling
out <- smote_oversampling(taxonomies_table, kmerMatrix, taxa)
kmerMatrix <- out[[1]]
taxonomies_table <- out[[2]]

# smote_rows 
smote_rows <- which(!startsWith(rownames(kmerMatrix), 'ID-'))
smote_data <- kmerMatrix[smote_rows, ]
smote_y <- taxonomies_table[smote_rows,]

# to be added in the training set
Xy_to_add <- cbind(smote_data,smote_y)

# delete smote rows
if (length(smote_rows) > 0){
  kmerMatrix <- kmerMatrix[-smote_rows,]
  taxonomies_table <- taxonomies_table[-smote_rows,]
}

print(length(smote_rows))

# Train test split
Xy <- cbind(kmerMatrix, taxonomies_table)

# find dataset characterisation metrics
binary_y <- y_binary_vector_encoding(taxonomies_table, 
                                     Xy[, c((ncol(Xy) - length(taxa) + 1):ncol(Xy))], 
                                     taxa)

binary_Xy <- cbind(Xy[, -c((ncol(Xy) - length(taxa) + 1):ncol(Xy))], binary_y)

#my_mldr <- mldr_from_dataframe(binary_Xy,  
#                               labelIndices = c((ncol(Xy) - length(taxa) + 1):ncol(binary_Xy)),
#                               name = "testMLDR")
#
#outputs <- data.table(measures = names(my_mldr$measures),
#                      values = unlist(my_mldr$measures))
#write.csv(outputs, 
#          paste0(output_folder, "/dataset characterisation metrics.csv"),
#          row.names = F)

rm(binary_Xy, binary_y)
rm(my_mldr, outputs)

# hyperparameter optimization
output_txt <- paste0(output_folder, "/grid-search-hp-opt-metrics.txt")
file.create(output_txt)
write("mtry ntrees accuracy micro-precision macro-precision micro-recall macro-recall micro-fmeasure macro-fmeasure hamming-loss subset-accuracy average_precision one_error coverage ranking_loss", output_txt, append = T)

class_opt <- floor(sqrt(ncol(kmerMatrix)))

mat <- rf_hyperparameters_optimization(Xy = Xy,
                                       Xy_to_add = Xy_to_add,
                                       num_of_experiments = num_of_experiments,
                                       mtry_list = mtry_list,
                                       ntrees_list = ntrees_list,
                                       taxa = taxa,
                                       taxonomies_table = taxonomies_table,
                                       output_txt = output_txt,
                                       user_threshold = predict_threshold)

#opt_results <- which(mat$metric == max(mat$metric))
#opt_results <- opt_results[1]
#mtry_opt <- mat$mtry[opt_results] 
#ntrees_opt <- mat$ntrees[opt_results] 

write.csv(mat,
          paste0(output_folder, "/hp grid search.csv"),
          row.names = FALSE)

#write.csv(mat[opt_results, ],
#            paste0(output_folder, "/optimal hyperparameters.csv"),
#            row.names = FALSE)


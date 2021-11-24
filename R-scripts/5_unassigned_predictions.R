# Clear
cat("\014")
rm(list = ls())

#libraries
library(data.table)
library(mlr)
library(mldr)
library(stringr)

# source
source('R-scripts/multilabel_functions.R')

# Set inputs -------------

# taxa
prediction_theshold <- 0.2

# input data
input_data <- read.csv('emp-data/Unassigned-sequences.csv') #[,-1]

# loading model 
my_model <- readRDS("Output/final-model/model.rds")

labels <- my_model$learner.model$yvar.names
y_test <- matrix(0,
                 nrow = nrow(input_data),
                 ncol = length(labels))

colnames(y_test) <- labels
y_test <- as.data.frame(y_test)

#load features
kmerMatrix <- read.csv("Output/preprocessing/kmerMatrix_fs.csv", row.names = 1)
kmers <- colnames(kmerMatrix)
rm(kmerMatrix)  

# data encoding -----------
X_test <- matrix(0L, 
                 nrow = nrow(input_data), 
                 ncol = length(kmers))

colnames(X_test) <- kmers
rownames(X_test) <- row.names(input_data)

for (one_kmer in kmers){
  X_test[, one_kmer] <- str_count(input_data$sequence, one_kmer)
}

X_test <- as.data.frame(X_test)
Xy_test <- cbind(X_test, y_test)

# Model predictions -----------

folder_model <- paste0("Output/unassigned-predictions") 
dir.create(folder_model)

# predict
predictions = predict(my_model, newdata = Xy_test)

# evaluation
# thr <- rep(prediction_theshold, length(labels))
# names(thr) <- labels
# predictions <- setThreshold(predictions, threshold = thr)


one_predictions <- keep_one_label(predictions$data)

fwrite(one_predictions$final_labels ,
       paste0(folder_model, "/unassigned-labels.csv"))


fwrite(one_predictions$y_pred ,
       paste0(folder_model, "/unassigned-predictions.csv"))



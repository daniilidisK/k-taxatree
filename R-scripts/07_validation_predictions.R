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
library(plyr)
library(dplyr)
library(parallel)
library(hash)
library(stringr)


# source
source('R-scripts/multilabel_functions.R')

# Prepare input to the model (train set) -------------

start_time <- Sys.time()
print(start_time)

# taxa
taxa <- c('kingdom', 'phylum', 'class', 'order')
prediction_theshold <- 0.2

my_model <- readRDS("Output/final-model/model.rds")
labels <- my_model$learner.model$yvar.names

# loading data matrix
# load features
kmerMatrix <- read.csv("Output/preprocessing/kmerMatrix_fs.csv", row.names = 1)
kmers <- colnames(kmerMatrix)
rm(kmerMatrix)  

# Prepare test set --------------
# loading testing data
taxonomies_table_testing <- read.csv('emp-data/emp-taxonomy-validation.csv')

# again, remove unassigned from the validation set
to_drop <- which(taxonomies_table_testing$kingdom == 'Unassigned')

# unassigned_seqs <- rbind(unassigned_seqs, 
#                          data.table(sequence = taxonomies_table_testing$sequence[to_drop]))
# 
# write.csv(unassigned_seqs, "Unassigned_sequences.csv", row.names = F)

if (length(to_drop) > 0){
  taxonomies_table_testing <- taxonomies_table_testing[-to_drop,]
}

testing_sequences <- taxonomies_table_testing$sequence

# data encoding -------------------
X_test <- matrix(0, 
                 nrow = nrow(taxonomies_table_testing), 
                 ncol = length(kmers))

colnames(X_test) <- kmers
rownames(X_test) <- taxonomies_table_testing$ID


for (one_kmer in kmers){
  X_test[, one_kmer] <- str_count(taxonomies_table_testing$sequence, one_kmer)
}

y_test <- taxonomies_table_testing[,taxa]

# encoding - taxonomy labels to binary
y_test_enc <- y_binary_vector_encoding(taxonomies_table_testing, y_test, taxa)

# colnames modification
# colnames need to be in a specified format (format of names of variables in R)
colnames(y_test_enc) <- gsub("-", "_", colnames(y_test_enc))      # replace '-' with '_'
colnames(y_test_enc) <- gsub('\\[','_lp_', colnames(y_test_enc))  # replace '[' with '_lp_' (left parenthesis)
colnames(y_test_enc) <- gsub('\\]','_rp_', colnames(y_test_enc))  # replace ']' with '_rp_' (right parenthesis)

y_test_enc <- y_test_enc[, labels]

y_test_enc <- as.data.frame(y_test_enc)
X_test <- as.data.frame(X_test)
Xy_test <- cbind(X_test, y_test_enc)


# Model predictions ----------------
folder_model <- paste0("Output/final-model") 
dir.create(folder_model)

# predict
predictions = predict(my_model, newdata = Xy_test)

# 
# thr <- rep(prediction_theshold, length(labels))
# names(thr) <- labels
# predictions <- setThreshold(predictions, threshold = thr)

#one_predictions <- keep_one_label(predictions$data)

#fwrite(one_predictions$final_labels ,
#       paste0(folder_model, "/labels.csv"))
#

#fwrite(one_predictions$y_pred ,
#       paste0(folder_model, "/predictions.csv"))

#y_pred <- one_predictions$y_pred

y_pred <- predictions$data

#y_pred <- y_pred[,which(startsWith(colnames(y_pred), 'prob.'))]
#colnames(y_pred) <- gsub("prob.","", colnames(y_pred))

# get predictions in [0,1] range
test_DT <- Xy_test  #data[test.set, ]

who <- which(stringr::str_detect(colnames(y_pred), "prob"))
prob_pred <- y_pred[, who]
colnames(prob_pred) <- stringr::str_remove_all(colnames(prob_pred), "prob.")

who <- which(stringr::str_detect(colnames(y_pred), "response"))
response_pred <- y_pred[, who]
colnames(response_pred) <- stringr::str_remove_all(colnames(response_pred), "response.")
response_pred[response_pred == TRUE] <- 1
response_pred[response_pred == FALSE] <- 0

response_pred <- sapply(response_pred, as.numeric)
# get all metrics 
return.table <- data.table(accuracy = 0,
                           micro_precision = 0,
                           macro_precision = 0,
                           micro_recall = 0,
                           macro_recall = 0,
                           micro_fmeasure = 0,
                           macro_fmeasure = 0) #,
                           #hamming_loss = 0,
                           #subset_accuracy = 0,
                           #average_precision = 0,
                           #one_error = 0,
                           #coverage = 0)

return.table$accuracy[1] <- accuracy(test_DT[, labels], 
                                     response_pred, 
                                     undefined_value = "ignore")

return.table$micro_precision[1] <- micro_precision(test_DT[, labels], 
                                                   response_pred, 
                                                   undefined_value = "ignore")

return.table$macro_precision[1] <- macro_precision(test_DT[, labels], 
                                                   response_pred, 
                                                   undefined_value = "ignore")

return.table$micro_recall[1] <- micro_recall(test_DT[, labels], 
                                             response_pred, 
                                             undefined_value = "ignore")

return.table$macro_recall[1] <- macro_recall(test_DT[, labels], 
                                             response_pred, 
                                             undefined_value = "ignore")

return.table$micro_fmeasure[1] <- micro_fmeasure(test_DT[, labels], 
                                                 response_pred, 
                                                 undefined_value = "ignore")

return.table$macro_fmeasure[1] <- macro_fmeasure(test_DT[, labels], 
                                                 response_pred, 
                                                 undefined_value = "ignore")

#return.table$hamming_loss[1] <- hamming_loss(test_DT[, labels], 
#                                             response_pred)
#
#return.table$subset_accuracy[1] <- subset_accuracy(test_DT[, labels], 
#                                                   response_pred)
#
#return.table$average_precision[1] <- average_precision(test_DT[, labels], 
#                                                       response_pred)
#
#return.table$one_error[1] <- one_error(test_DT[, labels], 
#                                       response_pred)
#
#return.table$coverage[1] <- coverage(test_DT[, labels], 
#                                     response_pred)
#
#return.table$ranking_loss[1] <- ranking_loss(test_DT[, labels], 
#                                             response_pred)

write.csv(return.table, 
          file = paste0(folder_model, '/final-metrics.csv'))



end_time <- Sys.time()

print(start_time - end_time)

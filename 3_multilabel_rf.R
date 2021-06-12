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

# models library
# library(C50)
# library(e1071)
#library(RWeka)

# source
source('multilabel_functions.R')

# taxa
taxa <- c('kingdom', 'phylum', 'class', 'order')

# put input data folder
input_folder <- "emp-data"

# initialization of output folders
output_folder <- 'multi-label-emp'
dir.create(output_folder)

for (tax in taxa){
  to_create <- paste0(output_folder, "/", tax)
  dir.create(to_create)
}

# taxonomies table
taxonomyFilepath <- paste0(input_folder, '/emp-taxonomy-train-test.csv')
taxonomies_table <- read.csv(taxonomyFilepath) #[,-1]
# taxonomies_table <- taxonomies_table[1:400, ]

# Creating data matrix
# nzv_corr file: data/dna_sequences_matrix_k_7.csv
# kmerAnalyzer file: data/kmerMatrix_nzv_corr.csv
kmerMatrix <- read.csv(paste0(input_folder, '/kmerMatrix.csv'),
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
kmerMatrix <- kmerMatrix[-smote_rows,]
taxonomies_table <- taxonomies_table[-smote_rows,]

# Train test split
Xy <- cbind(kmerMatrix, taxonomies_table)

# find dataset characterisation metrics
binary_y <- y_binary_vector_encoding(taxonomies_table, 
                                     Xy[, c((ncol(Xy) - length(taxa) + 1):ncol(Xy))], 
                                     taxa)

binary_Xy <- cbind(Xy[, -c((ncol(Xy) - length(taxa) + 1):ncol(Xy))], binary_y)

my_mldr <- mldr_from_dataframe(binary_Xy,  
                               labelIndices = c((ncol(Xy) - length(taxa) + 1):ncol(binary_Xy)),
                               name = "testMLDR")

outputs <- data.table(measures = names(my_mldr$measures),
                      values = unlist(my_mldr$measures))
write.csv(outputs, 
          paste0(output_folder, "/dataset characterisation metrics.csv"),
          row.names = F)

rm(binary_Xy, binary_y, my_mldr, outputs)

# hyperparameter optimization
# This might take some time. Instead, you can use the already optimized parameters 
# in lines 100, 101 (uncomment lines 94-95)

class_opt <- floor(sqrt(ncol(kmerMatrix)))
mtry_list <- seq(class_opt - 2, class_opt + 2, 1)
ntrees_list <- seq(700, 1000, 50)
num_of_experiments <- 5
mat <- rf_hyperparameters_optimization(Xy = Xy,
                                        Xy_to_add = Xy_to_add,
                                        num_of_experiments = num_of_experiments,
                                        mtry_list = mtry_list,
                                        ntrees_list = ntrees_list, 
                                        taxa = taxa,
                                        taxonomies_table = taxonomies_table)
 
opt_results <- which(mat$metric == max(mat$metric))
mtry_opt <- mat$mtry[opt_results] # for nzv_corr_boruta: 11, for kmerAnalyzer: 10
ntrees_opt <- mat$ntrees[opt_results] # for nzv_corr_boruta: 500, for kmerAnalyzer: 150

write.table(mat[opt_results, ], 
            paste0(output_folder, "/optimal hyperparameters.txt"),
            row.names = FALSE)

                  #anastasis      # using PCA_fs_0.9     #accuracy       #f1          #using PCA_fs_0.95    #f1    #using PCA_fs_1  #f1
#mtry_opt <- 21        #11                                   #10          13                                 14                      20
#ntrees_opt <- 100     #500                                 #300          500                                450                    400

# calculating metrics and feature importances
num_of_experiments = 5
out <- rf_average_metrics(Xy = Xy, 
                          Xy_to_add = Xy_to_add, 
                          num_of_experiments = num_of_experiments,
                          mtry_opt,
                          ntrees_opt,
                          taxa,
                          taxonomies_table,
                          output_folder)

# writing metrics to .csv file
eval_metrics <- out[[1]]
write.csv(eval_metrics, 
          file = paste0(output_folder,'/eval_metrics.csv'))

############################### TESTING ########################################

# train set
Xy_train <- rbind(Xy, Xy_to_add)

# Separating X_train and y_train
taxa_columns <- which(colnames(Xy_train) %in% taxa)
X_train <- Xy_train[,-taxa_columns]
y_train <- Xy_train[,taxa_columns]

# encoding - taxonomy labels to binary
y_train_enc <- y_binary_vector_encoding(taxonomies_table, y_train, taxa)

# merge again
Xy_train <- cbind(X_train, y_train_enc)

# Xy test
y_test <- matrix(FALSE, nrow = dim(X_test)[1], ncol = dim(y_train_enc)[2])
colnames(y_test) <- colnames(y_train_enc)
Xy_test <- cbind(X_test, y_test)

# colnames modification
# colnames need to be in a specified format (format of names of variables in R)
colnames(Xy_train) <- gsub("-", "_", colnames(Xy_train))      # replace '-' with '_'
colnames(Xy_train) <- gsub('\\[','_lp_', colnames(Xy_train))  # replace '[' with '_lp_' (left parenthesis)
colnames(Xy_train) <- gsub('\\]','_rp_', colnames(Xy_train))  # replace ']' with '_rp_' (right parenthesis)

colnames(Xy_test) <- gsub("-","_", colnames(Xy_test))
colnames(Xy_test) <- gsub('\\[','_lp_', colnames(Xy_test))
colnames(Xy_test) <- gsub('\\]','_rp_', colnames(Xy_test))

# finding again taxa columns in the new - encoded data
taxa_columns <- c( (dim(X_train)[2]+1) : dim(Xy_train)[2] )

# specify target columns
targets = colnames(Xy_train)[taxa_columns]

# merging data
data <- rbind(Xy_train, Xy_test)

n <- dim(X_train)[1]
train.set = seq(1, n, by = 1) # train rows
test.set = seq(n+1, dim(data)[1], by = 1) # test rows

# Multi label model
ml_task = makeMultilabelTask(data = data, target = targets)
classif.lrn <- makeLearner("multilabel.randomForestSRC", 
                           predict.type = "prob",
                           ntree = ntrees_opt, 
                           mtry = mtry_opt)

# training
model = train(classif.lrn, ml_task, subset = train.set)
saveRDS(model, 
        paste0(output_folder, "./final_model.rds"))

# predict
predictions = predict(model, 
                      task = ml_task, 
                      subset = test.set, 
                      type = 'prob')

thr <- rep(0.3, length(targets))
names(thr) <- targets
predictions <- setThreshold(predictions, threshold = thr)

y_pred <- predictions$data
y_pred <- y_pred[,which(startsWith(colnames(y_pred), 'response.'))]
colnames(y_pred) <- gsub("response.","", colnames(y_pred))

y_pred$label <- ""
y_pred <- as.data.table(y_pred)

for (i in 1:nrow(y_pred)) {
  
  temp <- melt.data.table(y_pred[i, ])[, value]
  who <- which(temp == TRUE)
  y_pred$label[i] <- paste(targets[who], collapse = ";")
}

# # evaluation
# # a <- performance(predictions, measures = list(multilabel.hamloss, multilabel.subset01, multilabel.ppv, multilabel.tpr, multilabel.f1, multilabel.acc))
# 
# # construct predictions
# y_pred_modified <- modify_predictions(y_pred, 
#                                       threshold = 0.1, 
#                                       taxa)
# 
# # create output table
# results_table <- y_pred_modified
# results_table$comparison <- results_table$`TRUE LABELS` == results_table$`PRED LABELS`
# results_table$sequence <- sequences[rownames(results_table),]

test_DT <- data[test.set, ]

who <- which(stringr::str_detect(colnames(predictions$data), "prob"))
prob_pred <- predictions$data[, who]
colnames(prob_pred) <- stringr::str_remove_all(colnames(prob_pred), "prob.")

who <- which(stringr::str_detect(colnames(predictions$data), "response"))
response_pred <- predictions$data[, who]
colnames(response_pred) <- stringr::str_remove_all(colnames(response_pred), "response.")
response_pred[response_pred == TRUE] <- 1
response_pred[response_pred == FALSE] <- 0

metrics_df <- data.table(metric = c("accuracy", "precision", "micro_precision",
                                    "macro_precision", "recall", "micro_recall",
                                    "macro_recall", "fmeasure", "micro_fmeasure",
                                    "macro_fmeasure", "hamming_loss", "subset_accuracy",
                                    "average_precision", "one_error", "coverage", 
                                    "ranking_loss"),
                         value = 0)
# Averaged metrics
metrics_df$value[1] <- accuracy(test_DT[, targets], response_pred, undefined_value = "ignore")
metrics_df$value[2] <- precision(test_DT[, targets], response_pred, undefined_value = "ignore")
metrics_df$value[3] <- micro_precision(test_DT[, targets], response_pred, undefined_value = "ignore")
metrics_df$value[4] <- macro_precision(test_DT[, targets], response_pred, undefined_value = "ignore")
metrics_df$value[5] <- recall(test_DT[,targets], response_pred, undefined_value = "ignore")
metrics_df$value[6] <- micro_recall(test_DT[, targets], response_pred, undefined_value = "ignore")
metrics_df$value[7] <- macro_recall(test_DT[, targets], response_pred, undefined_value = "ignore")
metrics_df$value[8] <- fmeasure(test_DT[, targets], response_pred, undefined_value = "ignore")
metrics_df$value[9] <- micro_fmeasure(test_DT[, targets], response_pred, undefined_value = "ignore")
metrics_df$value[10] <- macro_fmeasure(test_DT[, targets], response_pred, undefined_value = "ignore") 

# Basic metrics
metrics_df$value[11] <- hamming_loss(test_DT[, targets], response_pred)
metrics_df$value[12] <- subset_accuracy(test_DT[, targets], response_pred)

# Ranking based metrics
metrics_df$value[13] <- average_precision(test_DT[, targets], prob_pred)
metrics_df$value[14] <- one_error(test_DT[, targets], prob_pred)
metrics_df$value[15] <- coverage(test_DT[, targets], prob_pred)
metrics_df$value[16] <- ranking_loss(test_DT[, targets], prob_pred)

output_table <- data.table(predicted_label = y_pred[, c("label")],
                           true_label = "Unassigned")
# write to csv
write.csv(output_table,
          paste0(output_folder, '/test_labels.csv'),
          row.names = FALSE)

write.csv(metrics_df,
          paste0(output_folder, "/test_metrics.csv"),
          row.names = FALSE)

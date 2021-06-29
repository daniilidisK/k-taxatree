# Clear
cat("\014")
rm(list = ls())
# dev.off(dev.list()["RStudioGD"])

# libraries
library(UBL)
library(splitstackshape)
library(mlr)
library(stringr)
library(data.table)
library(caret)

# source
source('multilabel_functions.R')

####### FUNCTIONS ##################

# taxa
taxa <- c('kingdom', 'phylum', 'class', 'order')

# specify output folder
output_folder <- 'multi-label-emp'

# taxonomies table
taxonomyFilepath <- 'emp-data-loc/emp-taxonomy-train-test-loc.csv'
taxonomies_table <- read.csv(taxonomyFilepath) #[,-1]

# loading data matrix
kmerMatrix <- read.csv('emp-data-loc/kmerMatrix_fs.csv', row.names = 1)
rownames(kmerMatrix) <- taxonomies_table$ID

# Dropping totally Unassigned/ Unclassified data
out <- excluding_unassigned(kmerMatrix, taxonomies_table, taxa)
kmerMatrix <- out[[1]]
taxonomies_table <- out[[2]]

# keep only taxa columns
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

############################### TESTING ########################################
# loading testing data
taxonomies_table_testing <- read.csv('emp-data-loc/emp-taxonomy-validation-loc.csv')
to_drop <- which(taxonomies_table_testing$kingdom == 'Unassigned')
if (length(to_drop) > 0){
  taxonomies_table_testing <- taxonomies_table_testing[-to_drop,]
}
testing_sequences <- taxonomies_table_testing$sequence

# forming X_test
X_test <- matrix(0L, nrow = nrow(taxonomies_table_testing), ncol = ncol(kmerMatrix))
colnames(X_test) <- colnames(kmerMatrix)
rownames(X_test) <- paste(taxonomies_table_testing$ID,'-test', sep = '')

for (kmer in colnames(kmerMatrix)){
  
  X_test[,kmer] <- str_count(taxonomies_table_testing$sequence, paste0("(?=",kmer,")"))
}

y_test <- taxonomies_table_testing[,taxa]
y_test[which(y_test == 'Unassigned', arr.ind = TRUE)] <- ''
taxonomies_table_testing[which(taxonomies_table_testing == 'Unassigned', arr.ind = TRUE)] <- ''

# train set
Xy_train <- rbind(Xy, Xy_to_add)

# Separating X_train and y_train
taxa_columns <- which(colnames(Xy_train) %in% taxa)
X_train <- Xy_train[,-taxa_columns]
y_train <- Xy_train[,taxa_columns]

# encoding - taxonomy labels to binary
y_train_enc <- y_binary_vector_encoding(taxonomies_table, y_train, taxa)
y_test_enc <- y_binary_vector_encoding(taxonomies_table_testing, y_test, taxa)

y_train_test <- rbind.fill(as.data.frame(y_train_enc), as.data.frame(y_test_enc))
y_train_test[which(is.na(y_train_test), arr.ind = TRUE)] <- FALSE

# merge again
# Xy_train <- cbind(X_train, y_train_enc)

X_train_test <- rbind(X_train, X_test)

# merging data
data <- cbind(X_train_test, y_train_test)

# colnames modification
# colnames need to be in a specified format (format of names of variables in R)
colnames(data) <- gsub("-", "_", colnames(data))      # replace '-' with '_'
colnames(data) <- gsub('\\[','_lp_', colnames(data))  # replace '[' with '_lp_' (left parenthesis)
colnames(data) <- gsub('\\]','_rp_', colnames(data))  # replace ']' with '_rp_' (right parenthesis)

# finding again taxa columns in the new - encoded data
taxa_columns <- c( (dim(X_train)[2]+1) : dim(data)[2] )

# specify target columns
targets = colnames(data)[taxa_columns]


n <- dim(X_train)[1]
train.set = seq(1, n, by = 1) # train rows
test.set = seq(n+1, dim(data)[1], by = 1) # test rows

# Multi label model
# loading optimal parameters
params <- read.csv(paste0(output_folder, '/optimal hyperparameters.csv'), header = T)
ml_task = makeMultilabelTask(data = data, target = targets)
classif.lrn <- makeLearner("multilabel.randomForestSRC", predict.type = "prob", ntree = params$ntrees, mtry = params$mtry)

# training
model = train(classif.lrn, ml_task, subset = train.set)
saveRDS(model, "model.rds")

# predict
predictions = predict(model, task = ml_task, subset = test.set, type = 'prob')

#preds <- setThreshold(predictions, threshold = 0.3) #?????
y_pred <- predictions$data
y_pred <- y_pred[,which(startsWith(colnames(y_pred), 'prob.'))]
colnames(y_pred) <- gsub("prob.","", colnames(y_pred))

# evaluation
test_metrics <- performance(predictions, measures = list(multilabel.hamloss, multilabel.subset01, multilabel.ppv, multilabel.tpr, multilabel.f1, multilabel.acc))
write.csv(test_metrics, file = paste0(output_folder, '/validation_metrics_1.csv'))

# construct predictions
y_pred_modified <- modify_predictions(y_pred, threshold = 0, taxa)

# create output table
true_labels <- c()
for (i in 1:nrow(y_test)){
  to_append <- paste(y_test[i,][y_test[i,] != ""], collapse = ";")
  to_append <- gsub("-","_",to_append)
  true_labels <- rbind(true_labels, to_append)
}

#results_table$'TRUE LABELS' <- true_labels
check <- c(y_pred_modified$'PRED LABELS' == true_labels)

results_table <- data.frame(cbind(y_pred_modified$'PRED LABELS', true_labels, check))
colnames(results_table) <- c('Pred', 'True', 'Check')
results_table$Sequence <- testing_sequences
rownames(results_table) <- taxonomies_table_testing$ID

# write to csv
write.csv(results_table,paste0(output_folder, '/validation_labels.csv'))


############################## Metrics  ########################################
library(mldr)

y_pred$label <- y_pred_modified$`PRED LABELS`
y_pred <- as.data.table(y_pred)

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


write.csv(metrics_df,
          paste0(output_folder, "/validation_metrics_2.csv"),
          row.names = FALSE)
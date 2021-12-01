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

# taxa
taxa <- c('kingdom', 'phylum', 'class', 'order')
prediction_theshold <- 0.5

# taxonomies table
# taxonomyFilepath <- 'emp-data-loc/emp-taxonomy-train-test-loc.csv'
taxonomyFilepath <- 'emp-data/emp-taxonomy-train-test.csv'
taxonomies_table <- read.csv(taxonomyFilepath) #[,-1]

# loading data matrix
# kmerMatrix <- read.csv('emp-data-loc/kmerMatrix_fs.csv', row.names = 1)
kmerMatrix <- read.csv('Output/preprocessing/kmerMatrix_fs.csv', row.names = 1)
rownames(kmerMatrix) <- taxonomies_table$ID

# Dropping totally Unassigned/ Unclassified data
out <- excluding_unassigned(kmerMatrix, taxonomies_table, taxa)
kmerMatrix <- out[[1]]
taxonomies_table <- out[[2]]

unassigned_kmers <- out[[3]]
unassigned_seqs <- data.table(sequence = out[[4]])

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
smote_y <- taxonomies_table[smote_rows, ]

# to be added in the training set
Xy_to_add <- cbind(smote_data,smote_y)

# delete smote rows
kmerMatrix <- kmerMatrix[-smote_rows, ]
taxonomies_table <- taxonomies_table[-smote_rows, ]

# Train test split
Xy <- cbind(kmerMatrix, taxonomies_table)

# Prepare test set --------------
# loading testing data
# taxonomies_table_testing <- read.csv('emp-data-loc/emp-taxonomy-validation-loc.csv')
taxonomies_table_testing <- read.csv('emp-data/emp-taxonomy-validation.csv')

# again, remove unassigned from the validation set
to_drop <- which(taxonomies_table_testing$kingdom == 'Unassigned')

unassigned_seqs <- rbind(unassigned_seqs, 
                         data.table(sequence = taxonomies_table_testing$sequence[to_drop]))

write.csv(unassigned_seqs, "Unassigned_sequences.csv", row.names = F)

if (length(to_drop) > 0){
  taxonomies_table_testing <- taxonomies_table_testing[-to_drop,]
}

testing_sequences <- taxonomies_table_testing$sequence

# forming X_test
X_test <- matrix(0L, 
                 nrow = nrow(taxonomies_table_testing), 
                 ncol = ncol(kmerMatrix))

colnames(X_test) <- colnames(kmerMatrix)
rownames(X_test) <- paste(taxonomies_table_testing$ID,'-test',
                          sep = '')

for (kmer in colnames(kmerMatrix)){
  
  X_test[, kmer] <- str_count(taxonomies_table_testing$sequence, 
                              kmer)
}

# excluding sequences inconsistent with the train set labels
for (one_taxa in taxa) {
  
  who <- which(!(taxonomies_table_testing[, one_taxa] %in% taxonomies_table[, one_taxa]))
  taxonomies_table_testing[who, one_taxa] <- "Unassigned"
  
}

y_test <- taxonomies_table_testing[,taxa]

# remove again entirely "unassigned" sequences
who <- which(rowSums(y_test == "Unassigned") == 4)

if (length(who) > 0) {
  y_test <- y_test[-who, ]
  taxonomies_table_testing <- taxonomies_table_testing[-who, ]
}

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
# params <- read.csv('multi-label-emp/optimal hyperparameters.csv', header = T)
ml_task = makeMultilabelTask(data = data, 
                             target = targets)

classif.lrn <- makeLearner("multilabel.randomForestSRC", 
                           predict.type = "prob", 
                           ntree = 300, #300
                           mtry = 17)   #17

# training
model = train(classif.lrn, ml_task, subset = train.set)
folder_model <- paste0("Output/final-model") 
dir.create(folder_model)
saveRDS(model, paste0(folder_model, "/model.rds"))


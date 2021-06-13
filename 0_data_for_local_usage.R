# Clear
cat("\014")
rm(list = ls())

# libraries
library(splitstackshape)
library(data.table)
library(seqinr)

# load data
data <- read.csv('emp-data/emp-taxonomy-validation.csv')

# undersample
my_data <-data[which(data$order == 'o__Thiohalorhabdales' | data$order == 'o__Spirobacillales'| data$order == 'o__E2' | data$order == 'o__GIF10' | data$order == 'o__SSS58A'),]

# split data
taxa <- c('kingdom', 'phylum', 'class', 'order')
train_test <- stratified(my_data, 
                         taxa,
                         0.6, 
                         bothSets = TRUE,
                         keep.rownames = TRUE)

# train-test data
train_data <- as.data.table(train_test[[1]])
rownames(train_data) <- train_data$rn
train_data <- train_data[,-1]

# validation data
validation_data <- as.data.table(train_test[[2]])
rownames(validation_data) <- validation_data$rn
validation_data <- validation_data[,-1]

# storing to emp-data-loc

# train-test data
sequences <- train_data$sequence
seq_names <- train_data$ID
write.fasta(sequences = as.list(sequences),
            names = seq_names, 
            file.out = 'emp-data-loc/emp-train-test-loc.fasta',
            as.string = TRUE,
            nbchar = 50)

write.csv(train_data, 'emp-data-loc/emp-taxonomy-train-test-loc.csv', row.names = F)

# validation data
sequences <- validation_data$sequence
seq_names <- validation_data$ID
write.fasta(sequences = as.list(sequences),
            names = seq_names, 
            file.out = 'emp-data-loc/emp-validation-loc.fasta',
            as.string = TRUE,
            nbchar = 50)

write.csv(train_data, 'emp-data-loc/emp-taxonomy-validation-loc.csv', row.names = F)
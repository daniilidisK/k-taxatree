# Clear
cat("\014")
rm(list = ls())

# libraries
library(empdata)
library(stringr)
library(seqinr)
library(plyr)
library(splitstackshape)

# creating emp-data directory
dir.create('empdata')

# loading observation matrix
taxonomies <- observation_metadata_150bp
sequences <- biom_data_150bp$sequence
seq_names <- paste('ID-', c(1:length(sequences)), sep = '')

# writing fasta
write.fasta(as.list(sequences), names = seq_names, file.out = 'empdata/emp-total.fasta', as.string = T, nbchar = 50)

# colnames of taxonomies
taxonomies <- taxonomies[,-1]
names(taxonomies) <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
taxonomies[which(taxonomies == 'Unclassified', arr.ind = T)] <- 'Unassigned'

for (i in 1:ncol(taxonomies)){
  
  cur_tax <- taxonomies[,i]
  cur_tax[which(cur_tax == '' | str_length(cur_tax) == 3)] <- 'Unassigned'
  taxonomies[,i] <- cur_tax
  
}


#forming taxonomies
IDs_df <- data.frame(seq_names)
colnames(IDs_df) <- 'ID'
taxonomies <- cbind(IDs_df, taxonomies)

# adding sequences
taxonomies$sequence <- sequences
write.csv(taxonomies, 'empdata/emp-taxonomy-total.csv', row.names = F)

# creating train-test and validation set
taxa <- c('kingdom', 'phylum', 'class', 'order')
train_test <- stratified(taxonomies, taxa, 0.3, bothSets = TRUE, keep.rownames = TRUE)

# Xy_train to data frame, set row names, bla bla
train.test <- as.data.frame(train_test[[1]])
rownames(train.test) <- train.test$rn
train.test <- train.test[,-1]

# Xy_test to data frame, set row names, bla bla
validation <- as.data.frame(train_test[[2]])
rownames(validation) <- validation$rn
validation <- validation[,-1]

# saving
write.fasta(as.list(train.test$sequence), names = train.test$ID, file.out = 'empdata/emp-train-test.fasta', as.string = T, nbchar = 50)
write.fasta(as.list(validation$sequence), names = validation$ID, file.out = 'empdata/emp-validation.fasta', as.string = T, nbchar = 50)

write.csv(train.test, 'empdata/emp-taxonomy-train-test.csv', row.names = F)
write.csv(validation, 'empdata/emp-taxonomy-validation.csv', row.names = F)

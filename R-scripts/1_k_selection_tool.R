# Clear
cat("\014")
rm(list = ls())

# libraries
library(parallel)
library(data.table)
library(cluster)
library(Rfast)
library(plyr)

source("R-scripts/count_kmers_functions.R")
# INPUTS -----------------

# input files
taxa_file <- read.csv('emp-data-loc/emp-taxonomy-train-test-loc.csv')

# taxa
taxa <- c('kingdom', 'phylum', 'class', 'order')

# values of k
kvals <- c(4:12)

# output folder
dir.create("Output")
output_folder <- "Output/k-selection"
dir.create(output_folder)

#ANALYSIS ------------------------

# sequences, seqnames
sequences <- taxa_file$sequence
seq_names <- taxa_file$ID

# initialization
silh_scores_matrix <- matrix(0L, nrow = length(kvals), ncol = length(taxa))
rownames(silh_scores_matrix) <- kvals
colnames(silh_scores_matrix) <- taxa

output_txt <- paste0(output_folder, "/silhouette_scores.txt")
file.create(output_txt)
write(paste(taxa, collapse =","), output_txt, append = T)

# Analysis
for (k in kvals){
  
  print(paste0("Analysis for k = ", k))
  kmerMatrix <- count_kmers_in_file(sequences = sequences, seq_names = seq_names, k = k)
  print('Kmer matrix created')
  
  # calculate pairwise distances
  D <- Dist(kmerMatrix, method = 'minkowski', p = 1)
  
  print('Distance matrix calculated')
  rm(kmerMatrix)
  
  silh_row <- c()
  
  for (tax in taxa){
    
    print(paste0('Taxa: ', tax))
    
    targets <- taxa_file[,tax]
    to_drop <- which(targets == 'Unassigned')
    
    if(length(to_drop)>0){
      
      t <- factor(targets[-to_drop])
      s <- silhouette(as.numeric(t), D[-to_drop, -to_drop])
      silh_row <- c(silh_row, mean(s[,3]))
      
    } else {
      
      t <- factor(targets)
      s <- silhouette(as.numeric(t), D)
      silh_row <- c(silh_row, mean(s[,3]))
    
    }
  }
  
  silh_scores_matrix[as.character(k),] <- silh_row
  write(silh_row, output_txt, append = T)
  
  rm(silh_row, targets, D, t, s, to_drop)
  
}

results_txt <- paste0(output_folder, "/optimal k values.txt")

file.create(results_txt)

for (tax in taxa) {
  
  one_line <- paste0("Optimal k value for clustering of level ", tax, ":", 
                     row.names(silh_scores_matrix)[which.max(silh_scores_matrix[, tax])])
  
  write(one_line, results_txt , append = T)
  
}

library(ggplot2)

data_plot <- melt.data.table(as.data.table(silh_scores_matrix))
colnames(data_plot) <- c("taxa", "silhouette_score")
data_plot$k <- rep(kvals, 4)

png(paste0(output_folder, "/Silhouette scores.png"),  
    width = 7, height = 7, units = 'in', res = 300)

data_plot %>%
  ggplot(aes(x = k, y = silhouette_score, group = taxa, color = taxa)) +
  geom_line() +
  # scale_color(discrete = TRUE) +
  ggtitle("Silhouette scores") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 14, family = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.key = element_rect(fill = NA)) +
  scale_x_continuous(breaks = kvals)


dev.off()



library(taxize)
library(data.table)
library(annotate)

rm(list = ls())

test <- read.csv("../emp-data/Unassigned-sequences.csv")    

# part1
x = paste0(">", test$ID[1:25], "\n", test$sequence[1:25], "\n")
x = paste(x, collapse = "")

res <- blastSequences(x,                               #test$sequence[1],
		                  #hitListSize = 20,               #per sequence????
			                timeout = 360, 
			                as=c("data.frame"))

# part2
x = paste0(">", test$ID[26:50], "\n", test$sequence[26:50], "\n")
x = paste(x, collapse = "")

res <- rbind(res, blastSequences(x, timeout = 360, as=c("data.frame")))

colnames(res) <- str_replace_all(colnames(res), "-", "_")

# select best hit
# min e-value
# max bitscore
res <- as.data.table(res)

new_res <- res

new_res <- res[res[, .I[Hsp_evalue == min(Hsp_evalue)], 
                   by=Iteration_query_def]$V1]

new_res <- new_res[new_res[, .I[Hsp_bit_score == max(Hsp_bit_score)], 
                           by=Iteration_query_def]$V1]
#usethis::edit_r_environ() #ENTREZ_KEY='e84d6c2185b5d0b571782842c9b2bbdcba08'
#taxize::use_entrez()

titloi <-c("superkingdom", "phylum","class", "order", "family", 
           "genus", "species")

results_all <- data.table("id" = character(),
                          "superkingdom" = character(), 
                          "phylum" = character(),
                          "class" = character(), 
                          "order" = character(), 
                          "family" = character(), 
                          "genus" = character(),
                          "species" = character())

for (one_seq in unique(new_res$Iteration_query_def)){
  
  one_res <- new_res[which(new_res$Iteration_query_def == one_seq),]
  uid <- genbank2uid(one_res$Hit_accession)
  
  taxx <- classification(uid, db="ncbi")
  
  taxx.names <- rbindlist(taxx)
  taxx.names <- taxx.names[, -3]
  
  taxx.names <- taxx.names[which(taxx.names$rank %in% titloi), ]
  
  taxx.names <- taxx.names[, .(rank = rank, count = .N), by = name] 
  
  taxx.names <- unique(taxx.names)
  
  taxx.names$percentage <- taxx.names$count / nrow(one_res)
  
  taxx.names$result <- paste(taxx.names$name, taxx.names$percentage, sep = " ")
  
  taxx.names <- taxx.names[,.(final = paste(result, collapse = ";")),
                           by=rank]
  
  taxx.names <- taxx.names[order(match(rank, titloi))]
  
  one.row <- data.table(transpose(taxx.names))[-1, ]
  
  
  colnames(one.row) <- taxx.names$rank
  one.row$id <- one_res$Iteration_query_def[1]
  #one.row <- cbind(quey_id = , one.row)
  print(nrow(results_all) +1)
  results_all <- plyr::rbind.fill(results_all, one.row)
  
}



fwrite(results_all, file = "results_table.csv", sep = ",")

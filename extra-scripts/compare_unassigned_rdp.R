library(data.table)
library(stringr)
library(mldr)

blast_results <- fread("results_table.csv")

who_no <- which(blast_results$superkingdom == "")

rdp_results <- fread("RDP/unassigned -defaultRDPTaxonomy-predictFormat.txt")

# rdp_results <- rdp_results[, c("V1", "V6")]

identical(blast_results$query_id, rdp_results$`SEQ-ID`)

# rdp_boolean <- data.frame("Bacteria" = rdp_results$V6 == "Bacteria",
#                          "Archaea" = rdp_results$V6 == "Archaea")

rdp_boolean <- data.frame("Bacteria" = rdp_results$Bacteria,
                          "Archaea" = rdp_results$Archaea)
                          

blast_boolean <- data.frame("Bacteria" = str_detect(blast_results$superkingdom, "Bacteria"),
                            "Archaea" = str_detect(blast_results$superkingdom, "Archaea"))

blast_boolean[blast_boolean == TRUE] <- 1
blast_boolean[blast_boolean == FALSE] <- 0

rdp_boolean[rdp_boolean == TRUE] <- 1
rdp_boolean[rdp_boolean == FALSE] <- 0

blast_boolean <- blast_boolean[-who_no, ]
rdp_boolean <- rdp_boolean[-who_no, ]
# get all metrics 
return.table <- data.table(accuracy = 0,
                           micro_precision = 0,
                           macro_precision = 0,
                           micro_recall = 0,
                           macro_recall = 0,
                           micro_fmeasure = 0,
                           macro_fmeasure = 0)

return.table$accuracy[1] <- accuracy(blast_boolean, 
                                     rdp_boolean, 
                                     undefined_value = "ignore")

return.table$micro_precision[1] <- micro_precision(blast_boolean, 
                                                   rdp_boolean, 
                                                   undefined_value = "ignore")

return.table$macro_precision[1] <- macro_precision(blast_boolean, 
                                                   rdp_boolean, 
                                                   undefined_value = "ignore")

return.table$micro_recall[1] <- micro_recall(blast_boolean, 
                                             rdp_boolean, 
                                             undefined_value = "ignore")

return.table$macro_recall[1] <- macro_recall(blast_boolean, 
                                             rdp_boolean, 
                                             undefined_value = "ignore")

return.table$micro_fmeasure[1] <- micro_fmeasure(blast_boolean, 
                                                 rdp_boolean, 
                                                 undefined_value = "ignore")

return.table$macro_fmeasure[1] <- macro_fmeasure(blast_boolean, 
                                                 rdp_boolean, 
                                                 undefined_value = "ignore")


write.csv(return.table, file = "RDP/rdp-default-unassigned-metrics.csv")


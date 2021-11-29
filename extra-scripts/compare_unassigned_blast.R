library(data.table)
library(stringr)
library(mldr)

blast_results <- fread("../results_table.csv")

who_no <- which(blast_results$superkingdom == "")

rf_results <- fread("../unassigned-predictions/unassigned-labels.csv")

rf_results$max_label <- unlist(str_split(rf_results$superkingdom, " "))[c(T,F)]

# rf_boolean <- data.frame("Bacteria" = rf_results$prob.k__Bacteria == rf_results$max_label,
#                          "Archaea" = rf_results$prob.k__Archaea == rf_results$max_label)

rf_boolean <- data.frame("Bacteria" = rf_results$max_label == "Bacteria",
                          "Archaea" = rf_results$max_label == "Archaea")


blast_boolean <- data.frame("Bacteria" = str_detect(blast_results$superkingdom, "Bacteria"),
                            "Archaea" = str_detect(blast_results$superkingdom, "Archaea"))

blast_boolean[blast_boolean == TRUE] <- 1
blast_boolean[blast_boolean == FALSE] <- 0

rf_boolean[rf_boolean == TRUE] <- 1
rf_boolean[rf_boolean == FALSE] <- 0

blast_boolean <- blast_boolean[-who_no, ]
rf_boolean <- rf_boolean[-who_no, ]
# get all metrics 
return.table <- data.table(accuracy = 0,
                           micro_precision = 0,
                           macro_precision = 0,
                           micro_recall = 0,
                           macro_recall = 0,
                           micro_fmeasure = 0,
                           macro_fmeasure = 0)

return.table$accuracy[1] <- accuracy(blast_boolean, 
                                     rf_boolean, 
                                     undefined_value = "ignore")

return.table$micro_precision[1] <- micro_precision(blast_boolean, 
                                                   rf_boolean, 
                                                   undefined_value = "ignore")

return.table$macro_precision[1] <- macro_precision(blast_boolean, 
                                                   rf_boolean, 
                                                   undefined_value = "ignore")

return.table$micro_recall[1] <- micro_recall(blast_boolean, 
                                             rf_boolean, 
                                             undefined_value = "ignore")

return.table$macro_recall[1] <- macro_recall(blast_boolean, 
                                             rf_boolean, 
                                             undefined_value = "ignore")

return.table$micro_fmeasure[1] <- micro_fmeasure(blast_boolean, 
                                                 rf_boolean, 
                                                 undefined_value = "ignore")

return.table$macro_fmeasure[1] <- macro_fmeasure(blast_boolean, 
                                                 rf_boolean, 
                                                 undefined_value = "ignore")


write.csv(return.table, 
          file = "metrics.csv")


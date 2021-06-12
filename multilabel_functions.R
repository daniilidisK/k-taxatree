excluding_unassigned <- function(kmerMatrix, 
                                 taxonomies_table, 
                                 taxa){
  
  # keep only the first level (domain or kingdom)
  taxonomies <- taxonomies_table[,taxa[1]]
  
  # excluding totally unassigned rows
  to_drop <- which(taxonomies == 'Unassigned' | taxonomies == 'Unclassified')
  
  # # X test
  # X_test <- kmerMatrix[to_drop, ]
  # test_sequences <- taxonomies_table[to_drop,]$Sequence
  
  
  # drop
  if (length(to_drop)>1){
    taxonomies_table <- taxonomies_table[-to_drop,]
    kmerMatrix <- kmerMatrix[-to_drop, ]
  }
  
  return(list(kmerMatrix, taxonomies_table))
}


excluding_singletons <- function(kmerMatrix, 
                                 taxonomies_table, 
                                 taxa) {
  
  # reversing taxa table
  reverse_taxa <- rev(taxa)
  
  # initialization
  to_be_deleted <- c()
  
  for (i in 1:length(reverse_taxa)){
    
    # specify tax level
    curr_level <- reverse_taxa[i]
    
    # count classes and elements/ identifying singletons
    counts_table <- plyr::count(taxonomies_table[,curr_level])
    singletons <- counts_table[which(counts_table$freq == 1),]$x
    
    # deleting current level from singletons
    # example: if a singleton is: d_A; p_B; c_C; o_D we keep d_A; p_B; c_C;
    for (singleton in singletons){
      
      # find positions singletons in data
      index_of_singleton <- which(taxonomies_table[,curr_level] == singleton)
      total_row <- taxonomies_table[index_of_singleton,]
      
      # if current level is domain/kingdom, delete the entire row
      if (curr_level == taxa[1]){
        
        to_be_deleted <- c(to_be_deleted, index_of_singleton)
      
      }
      
      # else delete current level from singletons
      else{
        
        taxonomies_table[index_of_singleton,curr_level] <- ''
        
        for (j in (i+1):length(reverse_taxa)){
          
          previous_taxa <- reverse_taxa[j]
          counts_table_previous_taxa <- plyr::count(taxonomies_table[,previous_taxa])
          
          # identifying singletons of previous level
          index_of_row <- which(counts_table_previous_taxa$x == taxonomies_table[index_of_singleton, previous_taxa])
          
          # if there are no singletons, break
          if (counts_table_previous_taxa[index_of_row,]$freq > 1){
            break
          } else {
            
            taxonomies_table[index_of_singleton,previous_taxa] <- ''
            
            if (previous_taxa == taxa[1]){
              to_be_deleted <- c(to_be_deleted, index_of_singleton)
              break
            }
          }
        }
      }
    }
  }
  
  # exclude singleton rows
  if (length(to_be_deleted)>0){
    
    kmerMatrix <- kmerMatrix[-to_be_deleted,]
    taxonomies_table <- taxonomies_table[-to_be_deleted,]
  }
  
  return(list(kmerMatrix, 
              taxonomies_table))
}


smote_oversampling <- function(taxonomies_table, kmerMatrix, taxa){
  
  reverse_taxa <- rev(taxa)
  
  for (level in reverse_taxa){
    
    # specifying y matrix
    y_labeled <- taxonomies_table[,level]
    assigned <- which(y_labeled != '')
    y_labeled <- y_labeled[assigned]
    
    # specifying X matrix
    X <- kmerMatrix[assigned,]
    
    # # count classes -> identify which classes to oversample
    # count_classes <- plyr::count(y_labeled)
    # minority_classes_1 <- count_classes[which(count_classes$freq <= 5),]$x
    # minority_classes_2 <- count_classes[which(count_classes$freq > 5 & count_classes$freq <= 10),]$x
    # 
    # # Specifying smote coefficients
    # smote_over_coeffs <- hash()
    # 
    # if (level != taxa[1]) {
    # 
    #   for (min_class in minority_classes_1){
    #     smote_over_coeffs[min_class] <- 4
    #   }
    #   
    #   for (min_class in minority_classes_2){
    #     smote_over_coeffs[min_class] <- 2
    #   }
    # } else {
    #   
    #   for (min_class in minority_classes_1){
    #     smote_over_coeffs[min_class] <- 8
    #   }
    #   
    #   for (min_class in minority_classes_2){
    #     smote_over_coeffs[min_class] <- 4
    #   }
    # }
    # 
    
    # count classes -> identify which classes to oversample
    count_classes <- plyr::count(y_labeled)
    minority_classes_1 <- count_classes[which(count_classes$freq <= 10),]$x
    minority_classes_2 <- count_classes[which(count_classes$freq > 10 & count_classes$freq <= 50),]$x
    minority_classes_3 <- count_classes[which(count_classes$freq > 50 & count_classes$freq <= 100),]$x
    
    # Specifying smote coefficients
    smote_over_coeffs <- list()
    
    for (min_class in minority_classes_1){
      smote_over_coeffs[min_class] <- 10
    }
    
    for (min_class in minority_classes_2){
      smote_over_coeffs[min_class] <- 2
    }
    
    for (min_class in minority_classes_3){
      smote_over_coeffs[min_class] <- 1
    } 
    
    # apply smote algorithm
    Xy <- cbind(X, Y_ = as.factor(y_labeled))
    smote_over_coeffs <- as.list(smote_over_coeffs)
    newData <- SmoteClassif(Y_ ~ ., Xy, C.perc = smote_over_coeffs, k=1, dist = "Euclidean")
    
    # renaming smote samples
    smote_rows <- which(!startsWith(rownames(newData), 'ID-') & !startsWith(rownames(newData), 'smote') )
    rownames(newData)[smote_rows] <- paste('smote', level, c(1:length(smote_rows)), sep = '_' )
    
    # place generated smote data at the bottom of X matrix
    X <- newData[,-dim(Xy)[2]]
    smote_data <- X[smote_rows,]
    kmerMatrix <- rbind(kmerMatrix, smote_data)
    
    # place generated smote data in the bottom of y vector
    y_level <- newData[, 'Y_']
    y_level <- vapply(y_level, paste, collapse = ", ", character(1L))
    y_smote <- y_level[smote_rows]
    
    # assigning full path to generated smote data
    # step 1: initialization
    smote_y_table <- matrix('', ncol = length(taxa), nrow = length(y_smote))
    rownames(smote_y_table) <- rownames(smote_data)
    colnames(smote_y_table) <- rev(reverse_taxa)
    
    # step 2: to data frame
    smote_y_table <- as.data.frame(smote_y_table)
    
    # step 3: assign current level values
    smote_y_table[, level] <- y_smote
    
    # step 4: assign the rest of levels
    smote_labels <- unique(y_smote)
    for (label in smote_labels){
      level_index <- which(colnames(smote_y_table) == level)
      full_path <- taxonomies_table[which(taxonomies_table[,level] == label),][1,]
      smote_y_table[which(smote_y_table[,level] == label),1:level_index] <- full_path[,1:level_index]
    }
    
    taxonomies_table <- rbind(taxonomies_table, smote_y_table)
  }
  
  return(list(kmerMatrix, taxonomies_table))
}


y_binary_vector_encoding <- function(taxonomies_table, 
                                     y, 
                                     taxa){
  
  # all names
  taxa_names <- c(unique(taxonomies_table[,taxa[1]]), 
                             unique(taxonomies_table[,taxa[2]]),
                             unique(taxonomies_table[,taxa[3]]),
                             unique(taxonomies_table[,taxa[4]]))
  
  # exclude unassigned
  taxa_names <- taxa_names[taxa_names != ""]
  
  # initialization
  y_mat <- matrix(FALSE, nrow = nrow(y), ncol = length(taxa_names))
  colnames(y_mat) <- taxa_names
  rownames(y_mat) <- rownames(y)
  
  # binary encoding
  for(i in 1:nrow(y)){
    for (j in 1:ncol(y)){
      
      if (y[i,j] != ""){
        y_mat[i, y[i,j]] <- TRUE
      }
    }
  }
  
  return(y_mat)
}



rf_hyperparameters_optimization <- function(Xy, 
                                            Xy_to_add, 
                                            num_of_experiments, 
                                            mtry_list,
                                            ntrees_list, 
                                            taxa,
                                            taxonomies_table) {
  
  # binary encoding
  taxa_columns <- which(colnames(Xy) %in% taxa)
  y <- Xy[,taxa_columns]
  Xy <- Xy[,-taxa_columns]
  n_kmers <- ncol(Xy)
  
  # Separating X_test and y_test
  y_to_add <- Xy_to_add[,taxa_columns]
  Xy_to_add <- Xy_to_add[,-taxa_columns]
  
  # encoding - taxonomy labels to binary
  y_enc <- y_binary_vector_encoding(taxonomies_table, y, taxa)
  y_to_add_enc <- y_binary_vector_encoding(taxonomies_table, y_to_add, taxa)
  
  # merge again
  Xy <- cbind(Xy, y_enc)
  Xy_to_add <- cbind(Xy_to_add, y_to_add_enc)
  
  colnames(Xy) <- gsub("-", "_", colnames(Xy))      # replace '-' with '_'
  colnames(Xy) <- gsub('\\[','_lp_', colnames(Xy))  # replace '[' with '_lp_' (left parenthesis)
  colnames(Xy) <- gsub('\\]','_rp_', colnames(Xy))  # replace ']' with '_rp_' (right parenthesis)
  
  colnames(Xy_to_add) <- gsub("-","_", colnames(Xy_to_add))
  colnames(Xy_to_add) <- gsub('\\[','_lp_', colnames(Xy_to_add))
  colnames(Xy_to_add) <- gsub('\\]','_rp_', colnames(Xy_to_add))
  
  # prepare input vectors
  inputs <- data.table::data.table(mtry = numeric(),
                                   ntrees = numeric(),
                                   metric = numeric())
  # create the model 
  for (counter in 1:num_of_experiments) {
    
    new_table <- create_rf(counter,
                           mtry_list,
                           ntrees_list,
                           Xy, 
                           Xy_to_add,
                           n_kmers)
    
    inputs <- rbind(inputs, new_table)
  }
  
  # 
  # inputs$av_f1 <- mcmapply(create_rf, 
  #                          inputs$counter, 
  #                          MoreArgs = list(mtry_list = mtry_list,
  #                                          ntrees_list = ntrees_list,
  #                                          Xy = Xy, 
  #                                          Xy_to_add = Xy_to_add, 
  #                                          n_kmers = n_kmers),
  #                          mc.cores = 4) 
  
 # inputs <- inputs[, -3]
  inputs <- group_by(inputs,mtry, ntrees) %>%
    summarise(mtry, ntrees, metric = sum(metric)) %>%
    unique()
  
  inputs$metric <- inputs$metric / num_of_experiments
  
  return(inputs)
  
}


create_rf <- function(counter,
                      mtry_list,
                      ntrees_list,
                      Xy,
                      Xy_to_add,
                      n_kmers) {
  
  print(paste0("Iteration: ", counter))
  
  taxa <- colnames(Xy)[(n_kmers+1):ncol(Xy)]
  train_test <- stratified(Xy, 
                           taxa,
                           0.6, 
                           bothSets = TRUE,
                           keep.rownames = TRUE)
  
  # Xy_train to data frame, set row names, bla bla
  Xy_train <- as.data.table(train_test[[1]])
  rownames(Xy_train) <- Xy_train$rn
  Xy_train <- Xy_train[,-1]
  
  # add smote rows to training set
  Xy_train <- rbind(Xy_train, Xy_to_add)
  
  # Xy_test to data frame, set row names, bla bla
  Xy_test <- as.data.table(train_test[[2]])
  rownames(Xy_test) <- Xy_test$rn
  Xy_test <- Xy_test[,-1]
  
  # Separating X_train and y_train
  X_train <- Xy_train[, 1:n_kmers]
  y_train <- Xy_train[, (n_kmers+1):ncol(Xy_train)]
  
  # Separating X_test and y_test
  X_test <- Xy_test[, 1:n_kmers]
  y_test <- Xy_test[, (n_kmers+1):ncol(Xy_train)]
  
  # merge again
  Xy_train <- cbind(X_train, y_train)
  Xy_test <- cbind(X_test, y_test)
  
  # finding again taxa columns in the new - encoded data
  taxa_columns <- c((n_kmers+1):ncol(Xy_train))
  
  # specify target columns
  targets = colnames(Xy_train)[taxa_columns]
  
  # merging data
  data <- rbind(Xy_train, Xy_test)
  data <- as.data.frame(data)  
  train.set = seq(1, nrow(X_train)) # train rows
  test.set = seq(nrow(X_train)+1, nrow(data)) # test rows
  
  # Multi label model
  ml_task = makeMultilabelTask(data = data, target = targets)
  
  classif.lrn <- makeLearner("multilabel.randomForestSRC", 
                             predict.type = "prob", 
                             ntree = ntrees_list[1], 
                             mtry = mtry_list[1])
  
  # training
  model = train(classif.lrn, ml_task, subset = train.set)
  
  # # predict
  predictions = predict(model, 
                        task = ml_task, 
                        subset = test.set, 
                        type = 'prob')
  
  # get predictions in [0,1] range
  test_DT <- data[test.set, ]
  
  who <- which(stringr::str_detect(colnames(predictions$data), "prob"))
  prob_pred <- predictions$data[, who]
  colnames(prob_pred) <- stringr::str_remove_all(colnames(prob_pred), "prob.")
  
  who <- which(stringr::str_detect(colnames(predictions$data), "response"))
  response_pred <- predictions$data[, who]
  colnames(response_pred) <- stringr::str_remove_all(colnames(response_pred), "response.")
  response_pred[response_pred == TRUE] <- 1
  response_pred[response_pred == FALSE] <- 0
  
  # # Averaged metrics
  # accuracy(test_DT[, targets], response_pred, undefined_value = "ignore")
  # precision(test_DT[, targets], response_pred, undefined_value = "ignore")
  # micro_precision(test_DT[, targets], response_pred, undefined_value = "ignore")
  # macro_precision(test_DT[, targets], response_pred, undefined_value = "ignore")
  # recall(test_DT[,targets], response_pred, undefined_value = "ignore")
  # micro_recall(test_DT[, targets], response_pred, undefined_value = "ignore")
  # macro_recall(test_DT[, targets], response_pred, undefined_value = "ignore")
  # fmeasure(test_DT[, targets], response_pred, undefined_value = "ignore")
  # micro_fmeasure(test_DT[, targets], response_pred, undefined_value = "ignore")
  # macro_fmeasure(test_DT[, targets], response_pred, undefined_value = "ignore") #####
  # 
  # # Basic metrics
  # hamming_loss(test_DT[, targets], response_pred)
  # subset_accuracy(test_DT[, targets], response_pred)
  # 
  # # Ranking based metrics
  # average_precision(test_DT[, targets], prob_pred)
  # one_error(test_DT[, targets], prob_pred)
  # coverage(test_DT[, targets], prob_pred)
  # ranking_loss(test_DT[, targets], prob_pred)
  
  return.table <- data.table(mtry = rep(mtry_list, each = length(ntrees_list)),
                             ntrees = rep(ntrees_list, length(mtry_list)),
                             metric = 0)
  
  return.table$metric[1] <- macro_fmeasure(test_DT[, targets], 
                                           response_pred, 
                                           undefined_value = "ignore")
  
  for (j in 2:nrow(return.table)) {
    
    learner.new = setHyperPars(classif.lrn,
                               par.vals = list(ntree = return.table$ntrees[j], 
                                               mtry = return.table$mtry[j]))
    
    # training
    model.new = train(learner.new, ml_task, subset = train.set)
    
    # # predict
    predictions = predict(model.new, 
                          task = ml_task, 
                          subset = test.set, 
                          type = 'prob')
    
    who <- which(stringr::str_detect(colnames(predictions$data), "prob"))
    prob_pred <- predictions$data[, who]
    colnames(prob_pred) <- stringr::str_remove_all(colnames(prob_pred), "prob.")
    
    who <- which(stringr::str_detect(colnames(predictions$data), "response"))
    response_pred <- predictions$data[, who]
    colnames(response_pred) <- stringr::str_remove_all(colnames(response_pred), "response.")
    response_pred[response_pred == TRUE] <- 1
    response_pred[response_pred == FALSE] <- 0
    
    
    return.table$metric[j] <- macro_fmeasure(test_DT[, targets], 
                                             response_pred, 
                                             undefined_value = "ignore")
    
  }
  
  return(return.table) 
}


rf_average_metrics <- function(Xy, 
                               Xy_to_add,
                               num_of_experiments,
                               mtry_opt, 
                               ntrees_opt, 
                               taxa,
                               taxonomies_table,
                               output_folder) {
  
  # prepare the data
  # binary encoding
  taxa_columns <- which(colnames(Xy) %in% taxa)
  y <- Xy[,taxa_columns]
  Xy <- Xy[,-taxa_columns]
  n_kmers <- ncol(Xy)
  
  # Separating X_test and y_test
  y_to_add <- Xy_to_add[,taxa_columns]
  Xy_to_add <- Xy_to_add[,-taxa_columns]
  
  # encoding - taxonomy labels to binary
  y_enc <- y_binary_vector_encoding(taxonomies_table, y, taxa)
  y_to_add_enc <- y_binary_vector_encoding(taxonomies_table, y_to_add, taxa)
  
  # merge again
  Xy <- cbind(Xy, y_enc)
  Xy_to_add <- cbind(Xy_to_add, y_to_add_enc)
  
  colnames(Xy) <- gsub("-", "_", colnames(Xy))      # replace '-' with '_'
  colnames(Xy) <- gsub('\\[','_lp_', colnames(Xy))  # replace '[' with '_lp_' (left parenthesis)
  colnames(Xy) <- gsub('\\]','_rp_', colnames(Xy))  # replace ']' with '_rp_' (right parenthesis)
  
  colnames(Xy_to_add) <- gsub("-","_", colnames(Xy_to_add))
  colnames(Xy_to_add) <- gsub('\\[','_lp_', colnames(Xy_to_add))
  colnames(Xy_to_add) <- gsub('\\]','_rp_', colnames(Xy_to_add))
  
  metrics_df <- matrix(0, 
                       nrow = num_of_experiments,
                       ncol = 16)
  
  
  colnames(metrics_df) <- c("accuracy", "precision", "micro_precision",
                            "macro_precision", "recall", "micro_recall",
                            "macro_recall", "fmeasure", "micro_fmeasure",
                            "macro_fmeasure", "hamming_loss", "subset_accuracy",
                            "average_precision", "one_error", "coverage", 
                            "ranking_loss")
  importances <- hash()
  
  paths_prob <- list()
  
  for (counter in 1:num_of_experiments) {
    
    print(paste("Experiment ", counter, sep = '' ))
    
    taxa <- colnames(Xy)[(n_kmers+1):ncol(Xy)]
    
    train_test <- stratified(Xy, 
                             taxa,
                             0.6, 
                             bothSets = TRUE,
                             keep.rownames = TRUE)
    
    # Xy_train to data frame, set row names, bla bla
    Xy_train <- as.data.table(train_test[[1]])
    rownames(Xy_train) <- Xy_train$rn
    Xy_train <- Xy_train[,-1]
    
    # add smote rows to training set
    Xy_train <- rbind(Xy_train, Xy_to_add)
    
    # Xy_test to data frame, set row names, bla bla
    Xy_test <- as.data.table(train_test[[2]])
    rownames(Xy_test) <- Xy_test$rn
    Xy_test <- Xy_test[,-1]
    
    # Separating X_train and y_train
    X_train <- Xy_train[, 1:n_kmers]
    y_train <- Xy_train[, (n_kmers+1):ncol(Xy_train)]
    
    # Separating X_test and y_test
    X_test <- Xy_test[, 1:n_kmers]
    y_test <- Xy_test[, (n_kmers+1):ncol(Xy_train)]
    
    # merge again
    Xy_train <- cbind(X_train, y_train)
    Xy_test <- cbind(X_test, y_test)
    
    # finding again taxa columns in the new - encoded data
    taxa_columns <- c((n_kmers+1):ncol(Xy_train))
    
    # specify target columns
    targets = colnames(Xy_train)[taxa_columns]
    
    # merging data
    data <- rbind(Xy_train, Xy_test)
    data <- as.data.frame(data)  
    train.set = seq(1, nrow(X_train)) # train rows
    test.set = seq(nrow(X_train)+1, nrow(data)) # test rows
    
    # Multi label model
    ml_task = makeMultilabelTask(data = data, target = targets)
    classif.lrn <- makeLearner("multilabel.randomForestSRC", 
                               predict.type = "prob",
                               ntree = ntrees_opt, 
                               mtry = mtry_opt)
    
    # training
    model = train(classif.lrn, ml_task, subset = train.set)
    
    # feature importances
    if (counter == 1){
      vi <- randomForestSRC::vimp.rfsrc(model$learner.model)
      class_output <- vi$classOutput
      
      for(target in targets){
        temp <- class_output[target]
        importances[target] <- temp[target][[1]]$importance
      }
      
      importances <- as.list(importances)
      
      
    } else {
      
      vi <- randomForestSRC::vimp.rfsrc(model$learner.model)
      class_output <- vi$classOutput
      
      for(target in targets){
        
        temp <- class_output[target]
        which.name <- which(names(importances) == target)
        
        importances[[which.name]] <- importances[[which.name]] + temp[target][[1]]$importance
      }
    }
    
    # predict
    predictions = predict(model,
                          task = ml_task, 
                          subset = test.set,
                          type = 'prob')
    
    paths_prob[[counter]] <- predictions$data
    
    # get predictions in [0,1] range
    test_DT <- data[test.set, ]
    
    who <- which(stringr::str_detect(colnames(predictions$data), "prob"))
    prob_pred <- predictions$data[, who]
    colnames(prob_pred) <- stringr::str_remove_all(colnames(prob_pred), "prob.")
    
    who <- which(stringr::str_detect(colnames(predictions$data), "response"))
    response_pred <- predictions$data[, who]
    colnames(response_pred) <- stringr::str_remove_all(colnames(response_pred), "response.")
    response_pred[response_pred == TRUE] <- 1
    response_pred[response_pred == FALSE] <- 0
    
    # Averaged metrics
    metrics_df[counter, 1] <- accuracy(test_DT[, targets], response_pred, undefined_value = "ignore")
    metrics_df[counter, 2] <- precision(test_DT[, targets], response_pred, undefined_value = "ignore")
    metrics_df[counter, 3] <- micro_precision(test_DT[, targets], response_pred, undefined_value = "ignore")
    metrics_df[counter, 4] <- macro_precision(test_DT[, targets], response_pred, undefined_value = "ignore")
    metrics_df[counter, 5] <- recall(test_DT[,targets], response_pred, undefined_value = "ignore")
    metrics_df[counter, 6] <- micro_recall(test_DT[, targets], response_pred, undefined_value = "ignore")
    metrics_df[counter, 7] <- macro_recall(test_DT[, targets], response_pred, undefined_value = "ignore")
    metrics_df[counter, 8] <- fmeasure(test_DT[, targets], response_pred, undefined_value = "ignore")
    metrics_df[counter, 9] <- micro_fmeasure(test_DT[, targets], response_pred, undefined_value = "ignore")
    metrics_df[counter, 10] <- macro_fmeasure(test_DT[, targets], response_pred, undefined_value = "ignore") 
     
    # Basic metrics
    metrics_df[counter, 11] <- hamming_loss(test_DT[, targets], response_pred)
    metrics_df[counter, 12] <- subset_accuracy(test_DT[, targets], response_pred)
     
    # Ranking based metrics
    metrics_df[counter, 13] <- average_precision(test_DT[, targets], prob_pred)
    metrics_df[counter, 14] <- one_error(test_DT[, targets], prob_pred)
    metrics_df[counter, 15] <- coverage(test_DT[, targets], prob_pred)
    metrics_df[counter, 16] <- ranking_loss(test_DT[, targets], prob_pred)
    
    # new_mldr <- mldr_from_dataframe(test_DT[, targets], 
    #                                 labelIndices = c(1:length(targets)),
    #                                 name = "newMLDR")
    
    # cm <- multilabel_confusion_matrix(new_mldr, response_pred)
    # 
    # save_image(print(cm), 
    #            paste0(output_folder, "/cm_", counter, ".png"))
    
  }
  
  paths_prob <- rbindlist(paths_prob)
  
  
  write.table(paths_prob,  
            paste0(output_folder, "/paths_probs_table.csv"),
            row.names = FALSE)
  
  #rownames(metrics_df) <- c(1:num_of_experiments)
  
  for (target in targets){
    
    which.name <- which(names(importances) == target)
    importances[[which.name]] <- importances[[which.name]] / num_of_experiments
  }
  
  saving_feature_importances(importances, 
                             targets, 
                             output_folder)
  
  return(list(metrics_df, importances))
  
}


saving_feature_importances <- function(importances_list, 
                                       targets, 
                                       output_folder){
  
  for(target in targets){
    
    if (startsWith(target, 'k__')){
      temp_folder <- paste0(output_folder, "/kingdom")
    } else if (startsWith(target,'d__')){
      temp_folder <- paste0(output_folder, "/domain")
    } else if (startsWith(target, 'p__')){
      temp_folder <- paste0(output_folder, "/phylum")
    } else if (startsWith(target, 'c__')){
      temp_folder <- paste0(output_folder, "/class")
    } else if (startsWith(target, 'o__')){
      temp_folder <- paste0(output_folder, "/order")
    }
    
    which.name <- which(names(importances_list) == target)
    
    importances <- importances_list[[which.name]]
    
    write.csv(importances, 
              paste(temp_folder, '/', target,'.csv', sep = ''),
              row.names = TRUE)
  }
}


# there is also 'setThreshold' function from mlr
modify_predictions <- function(y_pred, threshold, taxa){
  
  # technically, tax_starts <- c('d__', 'p__', 'c__', 'o__')
  tax_starts <- paste(substring(taxa,1,1),'_', sep = '')
  
  # # colnames modification
  # y_test_modified <- data.frame(apply(y_test, 1, function(x) {
  #   u <- paste(x[x != ""], collapse = ";")
  #   u <- gsub("-","_",u) 
  #   u <- gsub('\\[','_lp_', u)
  #   u <- gsub('\\]','_rp_', u)
  #   u}))
  # 
  # # set colnames
  # colnames(y_test_modified) <- c('TRUE LABELS')
  # 
  # find predictions
  y_pred_modified <- data.frame(apply(y_pred, 1, function(x) {
    
    to_paste <- c()
    for (tax_start in tax_starts){
      to_append <- names(x)[which(startsWith(names(x), tax_start) & x>threshold)]
      #print(to_append)
      if (length(to_append) > 1){
        to_append <- to_append[which(x[to_append] == max(x[to_append]))]
        
        # if (length(to_append) > 1){
        #   print('hey op')
        # }
      }
      to_paste <- c(to_paste, to_append) 
    }
    
    paste(to_paste, collapse = ";")}))
  
  colnames(y_pred_modified) <- c('PRED LABELS')
  
  return(y_pred_modified)
  
}

## ------------------------------------------- ##
## FileName:  main_miRMarker
## Version: 0.1.0
## Description: This is a demo of miRMarker.
## Author:  DLUT AD&A Lab
## Created: 09-25 2023
## ------------------------------------------- ##


## Clear All Objects 
rm(list=ls(all.name=TRUE))
print('Clear All Objects')


## libraries 
setwd("E:\\Path\\Project")  
library(pROC)
library(mltools)
library(miRMarker)


## File Info 
file_name <- 'Files/GSE35834_01_mimat_unique.csv'

file_prefix <- miRMarker::GetFilePrefix(file_name)
start_t <- Sys.time()


## ------------------------------------------- ##

## Params (Cross Validation)
Times <- 10
Folds <- 10

## Params (Module Identification)
num_agent <- 5 
episodes_rl <- 2000
alpha_rl <- 0.2
gamma_rl <- 0.8
beta_net <- 0.5

## Params (Network)  
str_experi <- "ComDiffMinus"  # Cooperative Regulation Network

## Defined Classifier (FUN)
classifier <- miRMarker::SVMClassifier


## Expression Data
DataResults <- miRMarker::DataReader(file_name)
Datas <- DataResults[[1]]
Features <- DataResults[[2]]
Labels <- DataResults[[3]]

## Functional Similarity
file_disease_based <- "Files/miRNA_sim_matrix_disease_based.txt"
background_list_disease_based <- miRMarker::GetSimMatrixDatabase(file_disease_based)

## Results Dir
result_files_path <- miRMarker::MakeNewDir(getwd(), paste0('Results/', file_prefix))

## Initialization
list_matrix_test_auc <- miRMarker::ResultMatrixInit(1, Folds, Times)  # AUC 
list_matrix_test_sen <- miRMarker::ResultMatrixInit(1, Folds, Times)  # sensitivity
list_matrix_test_spe <- miRMarker::ResultMatrixInit(1, Folds, Times)  # specificity
list_matrix_test_mcc <- miRMarker::ResultMatrixInit(1, Folds, Times)  # MCC
matrix_module_weights <- matrix(0, nrow = Times*Folds, ncol = num_agent)  # module weights


## ------------------------------------------- ##

## Cross Validation
for(m in 1:Times)
{
  # Sample Division
  SampleFolds <- miRMarker::CrossValidationSampling(Labels, Folds, m)
  
  print("============================")
  
  # Cross Validation (single run)
  for(n in 1:Folds)
  {
    print("----------------------------")
    print(paste("Time", m, "Fold", n))
    
    # Training Samples, Test Samples
    test_Datas <- Datas[SampleFolds[[n]],]
    test_Labels <- Labels[SampleFolds[[n]]]
    train_Datas <- Datas[-SampleFolds[[n]],]
    train_Labels <- Labels[-SampleFolds[[n]]]
    
    # Data Scaling
    ScalingResults <- miRMarker::UVScaling(train_Datas, test_Datas)
    train_scaled <- ScalingResults[[1]]
    test_scaled <- ScalingResults[[2]]
    
    
    # Network Construction
    net_matrix <- miRMarker::NetworkExperimentalDatabaseWeighted(str_experi, beta_net, train_scaled, train_Labels, Features, 
                                                                 background_list_disease_based[["background_sim_matrix"]], 
                                                                 background_list_disease_based[["background_features"]])
    diag(net_matrix) <- 0.0
    
    
    # K Selection by PowerLaw Fitting
    k_value_list <- c(3, 5, 7, 9)
    net_matrix_final <- miRMarker::GraphMSTKNNDefined(net_matrix, k_value_list, "any")

    # Graph Info
    miRMarker::GraphInfo(net_matrix_final)
    
    
    # Agent Selection
    agent_list <- miRMarker::AgentSelectionByPvalueCentrality(net_matrix_final, train_scaled, train_Labels, "u_test", "closeness", num_agent, delta_value=0.5)
    
    
    # Module Identification by RL
    multi_module_list <- miRMarker::ModuleIdentifiedByRL(net_matrix_final, train_scaled, train_Labels, agent_list, classifier, episodes_rl, alpha_rl, gamma_rl)
    ModuleResults <- miRMarker::GetModuleDataFromModuleList(multi_module_list, net_matrix_final, train_scaled, test_scaled)
    
    
    # Prediction Probability (modules)
    num_modules <- length(ModuleResults[["node_index"]])
    prob_matrix_test <- matrix(0, nrow = length(test_Labels), ncol = num_modules)  # predicted probability
    for(k in 1:num_modules)
    {
      prob_matrix_test[, k] <- classifier(ModuleResults[["node_data_train"]][[k]], train_Labels, ModuleResults[["node_data_test"]][[k]])[[2]]
    }
    
    
    # Module Weights (weights)
    norm_type <- "softmax"
    module_weights <- miRMarker::ModuleWeightListCV(ModuleResults, train_Labels, classifier)  # weights
    module_index <- 1:num_modules
    module_index_sorted <- module_index[order(module_weights, decreasing = T)]  # sorted
    module_weights_sorted <- module_weights[order(module_weights, decreasing = T)]  # sorted
    matrix_module_weights[((m-1)*Folds+n), ] <- module_weights_sorted  # record sorted weights
    prob_matrix_test_temp <- matrix(prob_matrix_test[, module_index_sorted[1:num_agent]], ncol = num_agent)
    module_weights_temp <- module_weights_sorted[1:num_agent]
    
    
    # Weighted Averaging Predictions
    predict_prob_test_temp <- miRMarker::WeightedAverageModulePrediction(prob_matrix_test_temp, norm_type, module_weights_temp)  # predicted probability
    predict_labels_test_temp <- ifelse(predict_prob_test_temp >= 0.5, 1, 0)  # predicted labels
    write.csv(cbind(as.numeric(as.character(test_Labels)), predict_labels_test_temp, predict_prob_test_temp), 
              file=paste0(result_files_path, '/', file_prefix, '_', m, '_', n,'_a', num_agent,'.csv'), row.names=F)
    
    
    # Classification
    roc_model_test <- pROC::roc(test_Labels, predict_prob_test_temp, quiet=TRUE)  # ROC
    list_matrix_test_auc[[1]][n, m] <- roc_model_test$auc  # AUC
    print(roc_model_test$auc)
    
    sen_spe_results <- SenSpeForROC(roc_model_test)
    list_matrix_test_sen[[1]][n, m] <- sen_spe_results[[1]]  # Sensitivity
    list_matrix_test_spe[[1]][n, m] <- sen_spe_results[[2]]  # Specificity
    
    list_matrix_test_mcc[[1]][n, m] <- mltools::mcc(preds=predict_labels_test_temp, actuals=as.numeric(as.character(test_Labels)))  # MCC
  }
  
  print(Sys.time())
}


## ------------------------------------------- ##

## Save the Results
ResultMatrixSave(list_matrix_test_auc, paste0(file_prefix, '_cv_auc_', num_agent, '.csv'))
ResultMatrixSave(list_matrix_test_sen, paste0(file_prefix, '_cv_sen_', num_agent, '.csv'))
ResultMatrixSave(list_matrix_test_spe, paste0(file_prefix, '_cv_spe_', num_agent, '.csv'))
ResultMatrixSave(list_matrix_test_mcc, paste0(file_prefix, '_cv_mcc_', num_agent, '.csv'))
write.csv(matrix_module_weights, file = paste0(file_prefix, '_module_weights.csv'), row.names = F)


## Experiment Summary
exeTime <- difftime(Sys.time(), start_t, units="mins")
{
  print("--------------------------")
  print("miRMarker")
  print(file_name)
  print(paste('cv', Times, Folds))
  print("--------------------------")
  print(paste("Time:", as.numeric(exeTime), "mins"))
  print(Sys.time())
  print("--------------------------")
}
SaveResultSummary("miRMarker", file_prefix, Times, Folds, as.numeric(exeTime), as.character(Sys.time()))





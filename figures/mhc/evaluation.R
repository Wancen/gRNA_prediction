library(tidyverse)
library(data.table)

result_dir = "/proj/yunligrp/users/tianyou/gRNA/result/mhc"
xgboost_dir = "/proj/milovelab/mu/dukeproj/data/mhc/"
model_list = c("K562", "NPC", "iPSC")
CNN = list()
xgboost = list()
for (model in model_list){
  CNN_mat = matrix(NA, nrow = 5, ncol = 3)
  xgboost_mat = matrix(NA, nrow = 5, ncol = 3)
  for (fold in 1:5){
    for (test_cell in 1:3){
      CNN_result = fread(file.path(result_dir, tolower(model),
                                   paste0(tolower(model), "-binary-BCE-seqannot-fold", fold, 
                                          "-test-", tolower(model_list[test_cell]), ".csv")))
      if (model != model_list[test_cell]){
        xgboost_result = fread(file.path(xgboost_dir, tolower(model_list[test_cell]),"result",
                                         paste0(tolower(model_list[test_cell]),"-by-",tolower(model),
                                                "-pfdr0.05-pfdr0.2-binary-fold",fold,"-test.csv")),
                               col.names = c("xgboost_predict"))
      } else {
        xgboost_result = fread(file.path(xgboost_dir, tolower(model_list[test_cell]),"result",
                                         paste0(tolower(model_list[test_cell]),
                                                "-pfdr0.05-pfdr0.2-binary-fold",fold,"-test.csv")),
                               col.names = c("xgboost_predict"))
      }
      if (dim(xgboost_result)[1] != dim(CNN_result)[1]){
        warning(paste0("Model trained on ", model, " applied to ", model_list[celltype], " has inconsistent dimensions."))
      }
      result_merged = cbind(CNN_result, xgboost_result)
      CNN_mat[fold, test_cell] = pROC::auc(result_merged$true, result_merged$predict)
      xgboost_mat[fold, test_cell] = pROC::auc(result_merged$true, result_merged$xgboost_predict)
    }
  }
  CNN_mat = as.data.table(CNN_mat)
  colnames(CNN_mat) = model_list
  CNN[[model]] = CNN_mat
  xgboost_mat = as.data.table(xgboost_mat)
  colnames(xgboost_mat) = model_list
  xgboost[[model]] = xgboost_mat
}
CNN_combined = bind_rows(CNN, .id = "training_model")
write_csv(CNN_combined, file.path(result_dir, "mhc_auc_CNN_summary.csv"))

xgboost_combined = bind_rows(xgboost, .id = "training_model")
write_csv(xgboost_combined, file.path(result_dir, "mhc_auc_xgboost_summary.csv"))

combined = bind_rows(list(CNN=CNN_combined, XGBoost = xgboost_combined), .id="method")
write_csv(combined, file.path(result_dir, "mhc_auc_all_summary.csv"))

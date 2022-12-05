library(tidyverse)
library(data.table)

CNN_dir = "/proj/yunligrp/users/tianyou/gRNA/result/single-cell"
xgboost_dir = "/proj/milovelab/mu/dukeproj/data/scOct4th/result"
CNN_list = c("seq", "annot", "seqannot")
xgboost_list = c("-seq","-anno","")
savedir = "/proj/yunligrp/users/tianyou/gRNA/result/single-cell/comparison"
CNN_result_mat = matrix(NA, nrow = 5, ncol = length(CNN_list))
XGB_result_mat = matrix(NA, nrow = 5, ncol = length(CNN_list))

for (fold in 1:5){
  for (i in 1:length(CNN_list)){
    CNN_result = fread(file.path(CNN_dir, paste0("sc-binary-BCE-",
                                                 CNN_list[i],"-fold", fold, "-Dec04.csv")))
    xgb_result = fread(file.path(xgboost_dir, 
                                 paste0("sc_fdr_filtered-binary-",
                                        fold,"-test",xgboost_list[i],".csv")),
                       col.names = c("xgboost_predict"))
    if (dim(xgb_result)[1] != dim(CNN_result)[1]){
      warning(paste0("Model trained on fold ", fold, " with ", CNN_list[i], " has inconsistent dimensions."))
    }
    test_result = cbind(CNN_result, xgb_result)
    CNN_result_mat[fold, i] = pROC::auc(test_result$true, test_result$predict)
    XGB_result_mat[fold, i] = pROC::auc(test_result$true, test_result$xgboost_predict)
  }
}
CNN_result_mat = as.data.table(CNN_result_mat)
colnames(CNN_result_mat) = CNN_list
write_csv(CNN_result_mat, file.path(savedir, "sc_auc_CNN_summary.csv"))

XGB_result_mat = as.data.table(XGB_result_mat)
colnames(XGB_result_mat) = CNN_list
write_csv(XGB_result_mat, file.path(savedir, "sc_auc_xgboost_summary.csv"))

combined = bind_rows(list(CNN=CNN_result_mat, XGBoost = XGB_result_mat), .id="method")
write_csv(combined, file.path(savedir, "sc_auc_combined_summary.csv"))



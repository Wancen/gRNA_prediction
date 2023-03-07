library(tidyverse)
library(data.table)

CNN_dir = "/proj/yunligrp/users/tianyou/gRNA/result/binary_fivefold"
xgboost_dir = "/proj/milovelab/mu/dukeproj/data/dat_discovery/result"
CNN_list = c("seq", "seq-topannot", "seq-allannot", "annot")
savedir = "/proj/yunligrp/users/tianyou/gRNA/result/binary_fivefold/comparison"
xgboost_list = c("-seq","-top-seqanno","","-anno")
CNN_enh_result_mat = matrix(NA, nrow = 5, ncol = length(CNN_list))
XGB_enh_result_mat = matrix(NA, nrow = 5, ncol = length(CNN_list))

for (fold in 1:5){
  for (i in 1:length(CNN_list)){
    CNN_result = fread(file.path(CNN_dir, paste0("gRNA_binary-enh-BCE-",
                                                 CNN_list[i],"-fold", fold, "-Nov28.csv")))
    xgb_result = fread(file.path(xgboost_dir, 
                                 paste0("wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-",
                                        fold,"-test",xgboost_list[i],".csv")),
                       col.names = c("xgboost_predict"))
    if (dim(xgb_result)[1] != dim(CNN_result)[1]){
      warning(paste0("Model trained on fold ", fold, " with ", CNN_list[i], " has inconsistent dimensions."))
    }
    test_result = cbind(CNN_result, xgb_result)
    CNN_enh_result_mat[fold, i] = pROC::auc(test_result$true, test_result$predict)
    XGB_enh_result_mat[fold, i] = pROC::auc(test_result$true, test_result$xgboost_predict)
  }
}
CNN_enh_result_mat = as.data.table(CNN_enh_result_mat)
colnames(CNN_enh_result_mat) = CNN_list
write_csv(CNN_enh_result_mat, file.path(savedir, "enh_auc_CNN_summary.csv"))

XGB_enh_result_mat = as.data.table(XGB_enh_result_mat)
colnames(XGB_enh_result_mat) = CNN_list
write_csv(XGB_enh_result_mat, file.path(savedir, "enh_auc_xgboost_summary.csv"))

combined_enh = bind_rows(list(CNN=CNN_enh_result_mat, XGBoost = XGB_enh_result_mat), .id="method")
write_csv(combined_enh, file.path(savedir, "enh_auc_combined_summary.csv"))



##### For promoters ####
xgboost_list = c("-seq","-top-seqanno2","","-anno")
CNN_pro_result_mat = matrix(NA, nrow = 5, ncol = length(CNN_list))
XGB_pro_result_mat = matrix(NA, nrow = 5, ncol = length(CNN_list))

for (fold in 1:5){
  for (i in 1:length(CNN_list)){
    CNN_result = fread(file.path(CNN_dir, paste0("gRNA_binary-pro-BCE-",
                                                 CNN_list[i],"-fold", fold, "-Nov28.csv")))
    xgb_result = fread(file.path(xgboost_dir, 
                                 paste0("wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-",
                                        fold,"-test",xgboost_list[i],".csv")),
                       col.names = c("xgboost_predict"))
    if (dim(xgb_result)[1] != dim(CNN_result)[1]){
      warning(paste0("Model trained on fold ", fold, " with ", CNN_list[i], " has inconsistent dimensions."))
    }
    test_result = cbind(CNN_result, xgb_result)
    CNN_pro_result_mat[fold, i] = pROC::auc(test_result$true, test_result$predict)
    XGB_pro_result_mat[fold, i] = pROC::auc(test_result$true, test_result$xgboost_predict)
  }
}
CNN_pro_result_mat = as.data.table(CNN_pro_result_mat)
colnames(CNN_pro_result_mat) = CNN_list
write_csv(CNN_pro_result_mat, file.path(savedir, "pro_auc_CNN_summary.csv"))

XGB_pro_result_mat = as.data.table(XGB_pro_result_mat)
colnames(XGB_pro_result_mat) = CNN_list
write_csv(XGB_pro_result_mat, file.path(savedir, "pro_auc_xgboost_summary.csv"))

combined_pro = bind_rows(list(CNN=CNN_pro_result_mat, XGBoost = XGB_pro_result_mat), .id="method")
write_csv(combined_pro, file.path(savedir, "pro_auc_combined_summary.csv"))

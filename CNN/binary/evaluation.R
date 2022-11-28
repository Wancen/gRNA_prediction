library(tidyverse)
library(data.table)

result_dir = "/proj/yunligrp/users/tianyou/gRNA/result/binary_fivefold"
type_list = c("seq", "seq-topannot", "seq-allannot", "annot")
enh_result_mat = matrix(NA, nrow = 5, ncol = 4)
for (fold in 1:5){
  for (i in 1:4){
    test_result = fread(file.path(result_dir, paste0("gRNA_binary-enh-BCE-",
                                                     type_list[i],"-fold", fold, "-Nov28.csv")))
    enh_result_mat[fold, i] = pROC::auc(test_result$true, test_result$predict)
  }
}
enh_result_mat = as.data.table(enh_result_mat)
colnames(enh_result_mat) = type_list
write_csv(enh_result_mat, file.path(result_dir, "enh_auc_summary.csv"))

pro_result_mat = matrix(NA, nrow = 5, ncol = length(type_list))
for (fold in 1:5){
  for (i in 1:length(type_list)){
    test_result = fread(file.path(result_dir, paste0("gRNA_binary-pro-BCE-",
                                                     type_list[i],"-fold", fold, "-Nov28.csv")))
    pro_result_mat[fold, i] = pROC::auc(test_result$true, test_result$predict)
  }
}
pro_result_mat = as.data.table(pro_result_mat)
colnames(pro_result_mat) = type_list
write_csv(pro_result_mat, file.path(result_dir, "pro_auc_summary.csv"))


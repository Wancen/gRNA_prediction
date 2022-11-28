library(tidyverse)
library(data.table)

result_dir = "/proj/yunligrp/users/tianyou/gRNA/result_April_resplit/mhc"
celltype_list = c("k562", "npc", "ipsc")
result = list()
for (celltype in celltype_list){
  result_mat = matrix(NA, nrow = 5, ncol = 3)
  for (fold in 1:5){
    for (test_cell in 1:3){
      test_result = fread(file.path(result_dir, celltype, "binary", 
                                    paste0(celltype, "-binary-BCE-seqannot-fold", fold, "-test-", celltype_list[test_cell], ".csv")))
      result_mat[fold, test_cell] = pROC::auc(test_result$true, test_result$predict)
    }
  }
  result_mat = as.data.table(result_mat)
  colnames(result_mat) = celltype_list
  result[[celltype]] = result_mat
}


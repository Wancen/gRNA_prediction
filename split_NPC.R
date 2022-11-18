library(data.table)
library(tidyverse)

dir = "/proj/yunligrp/users/tianyou/gRNA/NPC"
dat = fread(file.path(dir, "NPC_onegene_full.csv"))
OGEE = fread("/proj/yunligrp/users/tianyou/gRNA/OGEE/NPC_OGEE.txt")
OGEE_downstream = fread("/proj/yunligrp/users/tianyou/gRNA/OGEE/NPC_OGEE_downstream.txt")
dat = dat %>% 
  left_join(OGEE %>% select(grna, distance, OGEE_prop_Essential), by = "grna") %>% 
  left_join(OGEE_downstream %>% select(grna, distance, OGEE_prop_Essential) %>% 
              rename(distance_downstream = distance,
                    OGEE_prop_Essential_downstream = OGEE_prop_Essential), by = "grna")
write_tsv(dat, file.path(dir, "split", "NPC_full.tsv"))

dat_sel = dat %>%
  filter(pval_fdr_corrected <= 0.05 | pval_fdr_corrected >= 0.2) %>%
  mutate(significant = ifelse(pval_fdr_corrected <= 0.05, 1, 0))
table(dat_sel$significant)
write_csv(dat_sel, file.path(dir, "split", "npc-rawp0.05.csv"))

## randomly sample dhs to split train and test data
## Binary classfication
set.seed(2050)
dhs <-dat_sel$dhs %>% unique() %>% sort()
sample = sample(1:length(dhs), length(dhs))
for (i in 1:5){
  s1 = ceiling(length(dhs)/5) * (i-1) + 1
  s2 = min(ceiling(length(dhs)/5) * i, length(dhs))
  dhs_fold = dhs[sample[s1:s2]]
  fold = dat_sel %>% filter(dhs %in% dhs_fold)
  print(paste0("Fold ", i, " # of gRNA: ", dim(fold)[1]))
  write_csv(fold, file.path(dir, "split", paste0("npc-pfdr0.05-pfdr0.2-binary-fold",i,".csv")))
}


library(tidyverse)
#library(DescTools)

######## For Promoters #######
dat = read_csv("wgCERES-gRNAs-k562-discovery-screen-pro_logFC0.1-test.csv")
colSums(is.na(dat))

dat_new = dat %>%
  select(-cancer_census_tissue_type, -cancer_census_role, -annotation_wg) %>%
  select(protospacer:chromEnd, gRNAid:DHS, baseMean:padj, everything())

dat_new = dat_new %>%
  mutate(ploidyZhou = ifelse(!is.na(ploidyZhou), ploidyZhou, 3),  ## imputation with mode
         distance = sign(distance) * log(abs(distance)+1),  ## log transformation
         medianRNAseqTPM = ifelse(!is.na(medianRNAseqTPM), medianRNAseqTPM, 0), ## log transformation for original values
         probIntolerantLoF = ifelse(!is.na(probIntolerantLoF), probIntolerantLoF, 0),
         numTKOHits_Hart = ifelse(!is.na(numTKOHits_Hart), numTKOHits_Hart, 0),
         OGEE_n_Essential = ifelse(!is.na(OGEE_n_Essential), OGEE_n_Essential, 0),
         OGEE_n_NonEssential = ifelse(!is.na(OGEE_n_NonEssential), OGEE_n_NonEssential, 0),
         OGEE_n = ifelse(!is.na(OGEE_n), OGEE_n, 0),
         OGEE_prop_Essential = ifelse(!is.na(OGEE_prop_Essential), OGEE_prop_Essential, 0),
         OGEE_prop_NonEssential = ifelse(!is.na(OGEE_prop_NonEssential), OGEE_prop_NonEssential, 0),
         vc_sqrt_sum = ifelse(!is.na(vc_sqrt_sum), vc_sqrt_sum, 0),
         DNAse_CPM_1Kb_new = ifelse(!is.na(DNAse_CPM_1Kb_new), DNAse_CPM_1Kb_new, 0),
         H3k27ac_CPM_1Kb_new = ifelse(!is.na(H3k27ac_CPM_1Kb_new), H3k27ac_CPM_1Kb_new, 0),
         ATAC_CPM_1Kb_new = ifelse(!is.na(ATAC_CPM_1Kb_new), ATAC_CPM_1Kb_new, 0),
         H3K4me3_CPM_1Kb_new = ifelse(!is.na(H3K4me3_CPM_1Kb_new), H3K4me3_CPM_1Kb_new, 0),
         TF_GATA2_CPM_1Kb_new = ifelse(!is.na(TF_GATA2_CPM_1Kb_new), TF_GATA2_CPM_1Kb_new, 0),
         TF_TAL1_CPM_1Kb_new = ifelse(!is.na(TF_TAL1_CPM_1Kb_new), TF_TAL1_CPM_1Kb_new, 0),
         TF_MYC_CPM_1Kb_new = ifelse(!is.na(TF_MYC_CPM_1Kb_new), TF_MYC_CPM_1Kb_new, 0))

colSums(is.na(dat_new))
write_csv(dat_new, "wgCERES-gRNAs-k562-discovery-screen-pro_logFC0.1-test-clean.csv")





# ####### For enhancers ########
dat_enh = read_csv("wgCERES-gRNAs-k562-discovery-screen-enh_logFC0.1-test.csv")
colSums(is.na(dat_enh))

dat_new_enh = dat_enh %>%
  select(-cancer_census_tissue_type, -cancer_census_role,-annotation_wg) %>%
  select(protospacer:chromEnd, gRNAid:DHS, baseMean:padj, everything())

dat_new_enh = dat_new_enh %>%
  mutate(ploidyZhou = ifelse(!is.na(ploidyZhou), ploidyZhou, 3),
         distance = sign(distance) * log(abs(distance)+1),  ## log transformation
         medianRNAseqTPM = ifelse(!is.na(medianRNAseqTPM), medianRNAseqTPM, 0), ## log transformation for original values
         probIntolerantLoF = ifelse(!is.na(probIntolerantLoF), probIntolerantLoF, 0),
         numTKOHits_Hart = ifelse(!is.na(numTKOHits_Hart), numTKOHits_Hart, 0),
         OGEE_n_Essential = ifelse(!is.na(OGEE_n_Essential), OGEE_n_Essential, 0),
         OGEE_n_NonEssential = ifelse(!is.na(OGEE_n_NonEssential), OGEE_n_NonEssential, 0),
         OGEE_n = ifelse(!is.na(OGEE_n), OGEE_n, 0),
         OGEE_prop_Essential = ifelse(!is.na(OGEE_prop_Essential), OGEE_prop_Essential, 0),
         OGEE_prop_NonEssential = ifelse(!is.na(OGEE_prop_NonEssential), OGEE_prop_NonEssential, 0),
         vc_sqrt_sum = ifelse(!is.na(vc_sqrt_sum), vc_sqrt_sum, 0),
         H3k27ac_CPM_1Kb_new = ifelse(!is.na(H3k27ac_CPM_1Kb_new), H3k27ac_CPM_1Kb_new, 0),
         DNAse_CPM_1Kb_new = ifelse(!is.na(DNAse_CPM_1Kb_new), DNAse_CPM_1Kb_new, 0),
         ATAC_CPM_1Kb_new = ifelse(!is.na(ATAC_CPM_1Kb_new), ATAC_CPM_1Kb_new, 0),
         H3K4me3_CPM_1Kb_new = ifelse(!is.na(H3K4me3_CPM_1Kb_new), H3K4me3_CPM_1Kb_new, 0),
         TF_GATA2_CPM_1Kb_new = ifelse(!is.na(TF_GATA2_CPM_1Kb_new), TF_GATA2_CPM_1Kb_new, 0),
         TF_TAL1_CPM_1Kb_new = ifelse(!is.na(TF_TAL1_CPM_1Kb_new), TF_TAL1_CPM_1Kb_new+1, 0),
         TF_MYC_CPM_1Kb_new = ifelse(!is.na(TF_MYC_CPM_1Kb_new), TF_MYC_CPM_1Kb_new, 0))

colSums(is.na(dat_new_enh))
write_csv(dat_new_enh, "wgCERES-gRNAs-k562-discovery-screen-enh_logFC0.1-test-clean.csv")




# 
# ### compare promoter vs enhancer distribution
# bind_rows(dat_new %>% select(protospacer, H3K27ac_CPM_per_1kbp), dat_new_enh %>% select(protospacer, H3K27ac_CPM_per_1kbp),
#           .id = "Pro/Enh") %>%
#   ggplot(aes(x = H3K27ac_CPM_per_1kbp, color = `Pro/Enh`)) +
#   geom_density()
# 
# 
### For promoters binary ####
dat = read_csv("wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-test.csv")
colSums(is.na(dat))

dat_new = dat %>%
  select(-cancer_census_tissue_type, -cancer_census_role,-annotation_wg) %>%
  select(protospacer:chromEnd, gRNAid:DHS, baseMean:padj, everything()) %>%
  select(-significance, everything())

dat_new = dat_new %>%
  mutate(ploidyZhou = ifelse(!is.na(ploidyZhou), ploidyZhou, 3),  ## imputation with mode
         distance = sign(distance) * log(abs(distance)+1),  ## log transformation
         medianRNAseqTPM = ifelse(!is.na(medianRNAseqTPM), medianRNAseqTPM, 0), ## log transformation for original values
         probIntolerantLoF = ifelse(!is.na(probIntolerantLoF), probIntolerantLoF, 0),
         numTKOHits_Hart = ifelse(!is.na(numTKOHits_Hart), numTKOHits_Hart, 0),
         OGEE_n_Essential = ifelse(!is.na(OGEE_n_Essential), OGEE_n_Essential, 0),
         OGEE_n_NonEssential = ifelse(!is.na(OGEE_n_NonEssential), OGEE_n_NonEssential, 0),
         OGEE_n = ifelse(!is.na(OGEE_n), OGEE_n, 0),
         OGEE_prop_Essential = ifelse(!is.na(OGEE_prop_Essential), OGEE_prop_Essential, 0),
         OGEE_prop_NonEssential = ifelse(!is.na(OGEE_prop_NonEssential), OGEE_prop_NonEssential, 0),
         vc_sqrt_sum = ifelse(!is.na(vc_sqrt_sum), vc_sqrt_sum, 0),
         DNAse_CPM_1Kb_new = ifelse(!is.na(DNAse_CPM_1Kb_new), DNAse_CPM_1Kb_new, 0),
         H3k27ac_CPM_1Kb_new = ifelse(!is.na(H3k27ac_CPM_1Kb_new), H3k27ac_CPM_1Kb_new, 0),
         ATAC_CPM_1Kb_new = ifelse(!is.na(ATAC_CPM_1Kb_new), ATAC_CPM_1Kb_new, 0),
         H3K4me3_CPM_1Kb_new = ifelse(!is.na(H3K4me3_CPM_1Kb_new), H3K4me3_CPM_1Kb_new, 0),
         TF_GATA2_CPM_1Kb_new = ifelse(!is.na(TF_GATA2_CPM_1Kb_new), TF_GATA2_CPM_1Kb_new, 0),
         TF_TAL1_CPM_1Kb_new = ifelse(!is.na(TF_TAL1_CPM_1Kb_new), TF_TAL1_CPM_1Kb_new, 0),
         TF_MYC_CPM_1Kb_new = ifelse(!is.na(TF_MYC_CPM_1Kb_new), TF_MYC_CPM_1Kb_new, 0))

colSums(is.na(dat_new))
write_csv(dat_new, "wgCERES-gRNAs-k562-discovery-screen-pro_baseMean125-binary-test-clean.csv")



### For enhancers binary ####
dat_enh = read_csv("wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-test.csv")
colSums(is.na(dat_enh))

dat_new_enh = dat_enh %>%
  select(-cancer_census_tissue_type, -cancer_census_role,-annotation_wg) %>%
  select(protospacer:chromEnd, gRNAid:DHS, baseMean:padj, everything()) %>%
  select(-significance, everything())

dat_new_enh = dat_new_enh %>%
  mutate(ploidyZhou = ifelse(!is.na(ploidyZhou), ploidyZhou, 3),
         distance = sign(distance) * log(abs(distance)+1),  ## log transformation
         medianRNAseqTPM = ifelse(!is.na(medianRNAseqTPM), medianRNAseqTPM, 0), ## log transformation for original values
         probIntolerantLoF = ifelse(!is.na(probIntolerantLoF), probIntolerantLoF, 0),
         numTKOHits_Hart = ifelse(!is.na(numTKOHits_Hart), numTKOHits_Hart, 0),
         OGEE_n_Essential = ifelse(!is.na(OGEE_n_Essential), OGEE_n_Essential, 0),
         OGEE_n_NonEssential = ifelse(!is.na(OGEE_n_NonEssential), OGEE_n_NonEssential, 0),
         OGEE_n = ifelse(!is.na(OGEE_n), OGEE_n, 0),
         OGEE_prop_Essential = ifelse(!is.na(OGEE_prop_Essential), OGEE_prop_Essential, 0),
         OGEE_prop_NonEssential = ifelse(!is.na(OGEE_prop_NonEssential), OGEE_prop_NonEssential, 0),
         vc_sqrt_sum = ifelse(!is.na(vc_sqrt_sum), vc_sqrt_sum, 0),
         H3k27ac_CPM_1Kb_new = ifelse(!is.na(H3k27ac_CPM_1Kb_new), H3k27ac_CPM_1Kb_new, 0),
         DNAse_CPM_1Kb_new = ifelse(!is.na(DNAse_CPM_1Kb_new), DNAse_CPM_1Kb_new, 0),
         ATAC_CPM_1Kb_new = ifelse(!is.na(ATAC_CPM_1Kb_new), ATAC_CPM_1Kb_new, 0),
         H3K4me3_CPM_1Kb_new = ifelse(!is.na(H3K4me3_CPM_1Kb_new), H3K4me3_CPM_1Kb_new, 0),
         TF_GATA2_CPM_1Kb_new = ifelse(!is.na(TF_GATA2_CPM_1Kb_new), TF_GATA2_CPM_1Kb_new, 0),
         TF_TAL1_CPM_1Kb_new = ifelse(!is.na(TF_TAL1_CPM_1Kb_new), TF_TAL1_CPM_1Kb_new+1, 0),
         TF_MYC_CPM_1Kb_new = ifelse(!is.na(TF_MYC_CPM_1Kb_new), TF_MYC_CPM_1Kb_new, 0))

colSums(is.na(dat_new_enh))
write_csv(dat_new_enh, "wgCERES-gRNAs-k562-discovery-screen-enh_baseMean125-binary-test-clean.csv")




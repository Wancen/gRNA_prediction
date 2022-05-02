library(tidyverse)
library(DescTools)

# ######## For Promoters #######
# dat = read_csv("wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05-padj0.2-train.csv")
# colSums(is.na(dat[,c(8:9,12:13,29:73)]))
# 
# dat_new = dat %>%
#   rename(strand = strand.y, H3K27ac_CPM_per_1kbp_new = H3k27ac_CPM_1Kb_new) %>%
#   select(protospacer:chromEnd, gRNAid:DHS, baseMean:padj, GCcount:Gcount, 
#          ploidyZhou:Gquad_n_overlap_other_strand,
#          distance:cancer_census_tier, vc_sqrt_sum:H3K27ac_CPM_per_1kbp, 
#          H3K27ac_CPM_per_1kbp_new:TF_MYC_CPM_1Kb_new,
#          deltagb:rna_dna, olapnumber:sumlog10fdr, strand, significant) %>%
#   mutate(strand = ifelse(strand == "+", 1, 
#                          ifelse(strand == "-", -1, NA)),
#          ploidyZhou = ifelse(!is.na(ploidyZhou), as.numeric(str_sub(ploidyZhou, 3, 3)), 3),
#          LossHetZhou = as.numeric(LossHetZhou),
#          SV_Zhou = as.numeric(SV_Zhou),
#          SNV_Zhou = as.numeric(SNV_Zhou),
#          Conserved = as.numeric(Conserved),
#          Gquad_same_strand = as.numeric(Gquad_same_strand),
#          Gquad_other_strand = as.numeric(Gquad_other_strand),
#          distance = sign(distance) * log(abs(distance)+1),  ## log transformation
#          medianRNAseqTPM = ifelse(!is.na(medianRNAseqTPM), medianRNAseqTPM, 0), ## log transformation for original values
#          probIntolerantLoF = ifelse(!is.na(probIntolerantLoF), probIntolerantLoF, 0),
#          numTKOHits_Hart = ifelse(!is.na(numTKOHits_Hart), numTKOHits_Hart, 0),
#          HartEssential = as.numeric(HartEssential),
#          OGEE_n_Essential = ifelse(!is.na(OGEE_n_Essential), OGEE_n_Essential, 0),
#          OGEE_n_NonEssential = ifelse(!is.na(OGEE_n_NonEssential), OGEE_n_NonEssential, 0),
#          OGEE_n = ifelse(!is.na(OGEE_n), OGEE_n, 0),
#          OGEE_prop_Essential = ifelse(!is.na(OGEE_prop_Essential), OGEE_prop_Essential, 0),
#          OGEE_prop_NonEssential = ifelse(!is.na(OGEE_prop_NonEssential), OGEE_prop_NonEssential, 0),
#          vc_sqrt_sum = ifelse(!is.na(vc_sqrt_sum), vc_sqrt_sum, 0),
#          DNase_CPM_per_1kbp = ifelse(!is.na(DNase_CPM_per_1kbp), DNase_CPM_per_1kbp, 0),
#          H3K27ac_CPM_per_1kbp = ifelse(!is.na(H3K27ac_CPM_per_1kbp), H3K27ac_CPM_per_1kbp, 0),
#          H3K27ac_CPM_per_1kbp_new = ifelse(!is.na(H3K27ac_CPM_per_1kbp_new), H3K27ac_CPM_per_1kbp_new, 0),
#          DNAse_CPM_1Kb_new = ifelse(!is.na(DNAse_CPM_1Kb_new), DNAse_CPM_1Kb_new, 0),
#          ATAC_CPM_1Kb_new = ifelse(!is.na(ATAC_CPM_1Kb_new), ATAC_CPM_1Kb_new, 0),
#          H3Kme4_CPM_1Kb_new = ifelse(!is.na(H3Kme4_CPM_1Kb_new), H3Kme4_CPM_1Kb_new, 0),
#          TF_GATA2_CPM_1Kb_new = ifelse(!is.na(TF_GATA2_CPM_1Kb_new), TF_GATA2_CPM_1Kb_new, 0),
#          TF_TAL1_CPM_1Kb_new = ifelse(!is.na(TF_TAL1_CPM_1Kb_new), TF_TAL1_CPM_1Kb_new, 0),
#          TF_MYC_CPM_1Kb_new = ifelse(!is.na(TF_MYC_CPM_1Kb_new), TF_MYC_CPM_1Kb_new, 0))
# 
# 
# write_csv(dat_new, "wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05-padj0.2-train-clean.csv")
# 
# 
# 
# 
# 
# ####### For enhancers ########
# dat_enh = read_csv("wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-padj0.2-test.csv")
# colSums(is.na(dat_enh[,c(8:13,29:73)]))
# 
# dat_new_enh = dat_enh %>%
#   rename(strand = strand.y, H3K27ac_CPM_per_1kbp_new = H3k27ac_CPM_1Kb_new) %>%
#   select(protospacer:chromEnd, gRNAid:DHS, baseMean:padj, GCcount:Gcount, 
#          ploidyZhou:Gquad_n_overlap_other_strand,
#          distance:cancer_census_tier, vc_sqrt_sum:H3K27ac_CPM_per_1kbp,
#          H3K27ac_CPM_per_1kbp_new:TF_MYC_CPM_1Kb_new,
#          deltagb:sumlog10fdr, strand, significant) %>%
#   mutate(strand = ifelse(strand == "+", 1, 
#                          ifelse(strand == "-", -1, NA)),
#          ploidyZhou = ifelse(!is.na(ploidyZhou), as.numeric(str_sub(ploidyZhou, 3, 3)), 3),
#          LossHetZhou = as.numeric(LossHetZhou),
#          SV_Zhou = as.numeric(SV_Zhou),
#          SNV_Zhou = as.numeric(SNV_Zhou),
#          Conserved = as.numeric(Conserved),
#          Gquad_same_strand = as.numeric(Gquad_same_strand),
#          Gquad_other_strand = as.numeric(Gquad_other_strand),
#          distance = sign(distance) * log(abs(distance)+1),  ## log transformation
#          medianRNAseqTPM = ifelse(!is.na(medianRNAseqTPM), log(medianRNAseqTPM+1), 0), ## log transformation for original values
#          probIntolerantLoF = ifelse(!is.na(probIntolerantLoF), probIntolerantLoF, 0),
#          numTKOHits_Hart = ifelse(!is.na(numTKOHits_Hart), numTKOHits_Hart, 0),
#          HartEssential = as.numeric(HartEssential),
#          OGEE_n_Essential = ifelse(!is.na(OGEE_n_Essential), OGEE_n_Essential, 0),
#          OGEE_n_NonEssential = ifelse(!is.na(OGEE_n_NonEssential), OGEE_n_NonEssential, 0),
#          OGEE_n = ifelse(!is.na(OGEE_n), OGEE_n, 0),
#          OGEE_prop_Essential = ifelse(!is.na(OGEE_prop_Essential), OGEE_prop_Essential, 0),
#          OGEE_prop_NonEssential = ifelse(!is.na(OGEE_prop_NonEssential), OGEE_prop_NonEssential, 0),
#          vc_sqrt_sum = ifelse(!is.na(vc_sqrt_sum), log(vc_sqrt_sum+1), 0),
#          DNase_CPM_per_1kbp = ifelse(!is.na(DNase_CPM_per_1kbp), log(DNase_CPM_per_1kbp+1), 0),
#          H3K27ac_CPM_per_1kbp = ifelse(!is.na(H3K27ac_CPM_per_1kbp), log(H3K27ac_CPM_per_1kbp+1), 0),
#          H3K27ac_CPM_per_1kbp_new = ifelse(!is.na(H3K27ac_CPM_per_1kbp_new), log(H3K27ac_CPM_per_1kbp_new+1), 0),
#          DNAse_CPM_1Kb_new = ifelse(!is.na(DNAse_CPM_1Kb_new), log(DNAse_CPM_1Kb_new+1), 0),
#          ATAC_CPM_1Kb_new = ifelse(!is.na(ATAC_CPM_1Kb_new), log(ATAC_CPM_1Kb_new+1), 0),
#          H3Kme4_CPM_1Kb_new = ifelse(!is.na(H3Kme4_CPM_1Kb_new), log(1+log(H3Kme4_CPM_1Kb_new+1)), 0),
#          TF_GATA2_CPM_1Kb_new = ifelse(!is.na(TF_GATA2_CPM_1Kb_new), log(TF_GATA2_CPM_1Kb_new+1), 0),
#          TF_TAL1_CPM_1Kb_new = ifelse(!is.na(TF_TAL1_CPM_1Kb_new), log(1+log(TF_TAL1_CPM_1Kb_new+1)), 0),
#          TF_MYC_CPM_1Kb_new = ifelse(!is.na(TF_MYC_CPM_1Kb_new), log(TF_MYC_CPM_1Kb_new+1), 0),
#          olapnumber = log(olapnumber + 1),
#          sumlog10fdr = log(1 + sumlog10fdr),
#          promlog10fdr = log(1 + promlog10fdr))
# 
# write_csv(dat_new_enh, "wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-padj0.2-test-clean.csv")
# 
# 
# ### compare promoter vs enhancer distribution
# bind_rows(dat_new %>% select(protospacer, H3K27ac_CPM_per_1kbp), dat_new_enh %>% select(protospacer, H3K27ac_CPM_per_1kbp),
#           .id = "Pro/Enh") %>%
#   ggplot(aes(x = H3K27ac_CPM_per_1kbp, color = `Pro/Enh`)) +
#   geom_density()
# 
# 
### For promoters binary ####
dat = read_csv("wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05-binary-train.csv")
colSums(is.na(dat[,c(8:9,12:13,29:73)]))

dat_new = dat %>%
  rename(strand = strand.y, H3K27ac_CPM_per_1kbp_new = H3k27ac_CPM_1Kb_new) %>%
  select(protospacer:chromEnd, gRNAid:DHS, baseMean:padj, GCcount:Gcount, 
         ploidyZhou:Gquad_n_overlap_other_strand,
         distance:cancer_census_tier, vc_sqrt_sum:H3K27ac_CPM_per_1kbp, 
         H3K27ac_CPM_per_1kbp_new:TF_MYC_CPM_1Kb_new,
         deltagb:rna_dna, olapnumber:sumlog10fdr, strand, significant)

for (i in c(29,39:48,51:52)){
  print(colnames(dat_new)[i])
  dat_new[,i] = Winsorize(dat_new[,i], probs = c(0,0.95), na.rm = T)
}

dat_new = dat_new %>%
  mutate(strand = ifelse(strand == "+", 1, 
                         ifelse(strand == "-", -1, NA)),
         ploidyZhou = ifelse(!is.na(ploidyZhou), as.numeric(str_sub(ploidyZhou, 3, 3)), 3),  ## imputation with mode
         LossHetZhou = as.numeric(LossHetZhou),
         SV_Zhou = as.numeric(SV_Zhou),
         SNV_Zhou = as.numeric(SNV_Zhou),
         Conserved = as.numeric(Conserved),
         Gquad_same_strand = as.numeric(Gquad_same_strand),
         Gquad_other_strand = as.numeric(Gquad_other_strand),
         distance = sign(distance) * log(abs(distance)+1),  ## log transformation
         medianRNAseqTPM = ifelse(!is.na(medianRNAseqTPM), medianRNAseqTPM, 0), ## log transformation for original values
         probIntolerantLoF = ifelse(!is.na(probIntolerantLoF), probIntolerantLoF, 0),
         numTKOHits_Hart = ifelse(!is.na(numTKOHits_Hart), numTKOHits_Hart, 0),
         HartEssential = as.numeric(HartEssential),
         OGEE_n_Essential = ifelse(!is.na(OGEE_n_Essential), OGEE_n_Essential, 0),
         OGEE_n_NonEssential = ifelse(!is.na(OGEE_n_NonEssential), OGEE_n_NonEssential, 0),
         OGEE_n = ifelse(!is.na(OGEE_n), OGEE_n, 0),
         OGEE_prop_Essential = ifelse(!is.na(OGEE_prop_Essential), OGEE_prop_Essential, 0),
         OGEE_prop_NonEssential = ifelse(!is.na(OGEE_prop_NonEssential), OGEE_prop_NonEssential, 0),
         vc_sqrt_sum = ifelse(!is.na(vc_sqrt_sum), vc_sqrt_sum, 0),
         DNase_CPM_per_1kbp = ifelse(!is.na(DNase_CPM_per_1kbp), DNase_CPM_per_1kbp, 0),
         H3K27ac_CPM_per_1kbp = ifelse(!is.na(H3K27ac_CPM_per_1kbp), H3K27ac_CPM_per_1kbp, 0),
         H3K27ac_CPM_per_1kbp_new = ifelse(!is.na(H3K27ac_CPM_per_1kbp_new), H3K27ac_CPM_per_1kbp_new, 0),
         DNAse_CPM_1Kb_new = ifelse(!is.na(DNAse_CPM_1Kb_new), DNAse_CPM_1Kb_new, 0),
         ATAC_CPM_1Kb_new = ifelse(!is.na(ATAC_CPM_1Kb_new), ATAC_CPM_1Kb_new, 0),
         H3Kme4_CPM_1Kb_new = ifelse(!is.na(H3Kme4_CPM_1Kb_new), H3Kme4_CPM_1Kb_new, 0),
         TF_GATA2_CPM_1Kb_new = ifelse(!is.na(TF_GATA2_CPM_1Kb_new), TF_GATA2_CPM_1Kb_new, 0),
         TF_TAL1_CPM_1Kb_new = ifelse(!is.na(TF_TAL1_CPM_1Kb_new), TF_TAL1_CPM_1Kb_new, 0),
         TF_MYC_CPM_1Kb_new = ifelse(!is.na(TF_MYC_CPM_1Kb_new), TF_MYC_CPM_1Kb_new, 0))

write_csv(dat_new, "wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05-binary-train-clean-winsorize.csv")



### For enhancers binary ####
dat_enh = read_csv("wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-binary-test.csv")
colSums(is.na(dat_enh[,c(8:13,29:73)]))

dat_new_enh = dat_enh %>%
  rename(strand = strand.y, H3K27ac_CPM_per_1kbp_new = H3k27ac_CPM_1Kb_new) %>%
  select(protospacer:chromEnd, gRNAid:DHS, baseMean:padj, GCcount:Gcount, 
         ploidyZhou:Gquad_n_overlap_other_strand,
         distance:cancer_census_tier, vc_sqrt_sum:H3K27ac_CPM_per_1kbp,
         H3K27ac_CPM_per_1kbp_new:TF_MYC_CPM_1Kb_new,
         deltagb:sumlog10fdr, strand, significant) 

for (i in c(29,39:48,51:54)){
  print(colnames(dat_new_enh)[i])
  dat_new_enh[,i] = Winsorize(dat_new_enh[,i], probs = c(0,0.95), na.rm = T)
}

dat_new_enh = dat_new_enh %>%
  mutate(strand = ifelse(strand == "+", 1, 
                         ifelse(strand == "-", -1, NA)),
         ploidyZhou = ifelse(!is.na(ploidyZhou), as.numeric(str_sub(ploidyZhou, 3, 3)), 3),
         LossHetZhou = as.numeric(LossHetZhou),
         SV_Zhou = as.numeric(SV_Zhou),
         SNV_Zhou = as.numeric(SNV_Zhou),
         Conserved = as.numeric(Conserved),
         Gquad_same_strand = as.numeric(Gquad_same_strand),
         Gquad_other_strand = as.numeric(Gquad_other_strand),
         distance = sign(distance) * log(abs(distance)+1),  ## log transformation
         medianRNAseqTPM = ifelse(!is.na(medianRNAseqTPM), medianRNAseqTPM, 0), ## log transformation for original values
         probIntolerantLoF = ifelse(!is.na(probIntolerantLoF), probIntolerantLoF, 0),
         numTKOHits_Hart = ifelse(!is.na(numTKOHits_Hart), numTKOHits_Hart, 0),
         HartEssential = as.numeric(HartEssential),
         OGEE_n_Essential = ifelse(!is.na(OGEE_n_Essential), OGEE_n_Essential, 0),
         OGEE_n_NonEssential = ifelse(!is.na(OGEE_n_NonEssential), OGEE_n_NonEssential, 0),
         OGEE_n = ifelse(!is.na(OGEE_n), OGEE_n, 0),
         OGEE_prop_Essential = ifelse(!is.na(OGEE_prop_Essential), OGEE_prop_Essential, 0),
         OGEE_prop_NonEssential = ifelse(!is.na(OGEE_prop_NonEssential), OGEE_prop_NonEssential, 0),
         vc_sqrt_sum = ifelse(!is.na(vc_sqrt_sum), vc_sqrt_sum, 0),
         DNase_CPM_per_1kbp = ifelse(!is.na(DNase_CPM_per_1kbp), DNase_CPM_per_1kbp, 0),
         H3K27ac_CPM_per_1kbp = ifelse(!is.na(H3K27ac_CPM_per_1kbp), H3K27ac_CPM_per_1kbp, 0),
         H3K27ac_CPM_per_1kbp_new = ifelse(!is.na(H3K27ac_CPM_per_1kbp_new), H3K27ac_CPM_per_1kbp_new, 0),
         DNAse_CPM_1Kb_new = ifelse(!is.na(DNAse_CPM_1Kb_new), DNAse_CPM_1Kb_new, 0),
         ATAC_CPM_1Kb_new = ifelse(!is.na(ATAC_CPM_1Kb_new), ATAC_CPM_1Kb_new, 0),
         H3Kme4_CPM_1Kb_new = ifelse(!is.na(H3Kme4_CPM_1Kb_new), H3Kme4_CPM_1Kb_new, 0),
         TF_GATA2_CPM_1Kb_new = ifelse(!is.na(TF_GATA2_CPM_1Kb_new), TF_GATA2_CPM_1Kb_new, 0),
         TF_TAL1_CPM_1Kb_new = ifelse(!is.na(TF_TAL1_CPM_1Kb_new), TF_TAL1_CPM_1Kb_new+1, 0),
         TF_MYC_CPM_1Kb_new = ifelse(!is.na(TF_MYC_CPM_1Kb_new), TF_MYC_CPM_1Kb_new, 0))

write_csv(dat_new_enh, "wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-binary-test-clean-winsorize.csv")




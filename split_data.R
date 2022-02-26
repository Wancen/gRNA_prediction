library(tidyverse)
require(rtracklayer)
# extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
#                           qValue = "numeric", peak = "integer")
# gr_narrowPeak <- import("C:/Users/wancenmu/OneDrive - University of North Carolina at Chapel Hill/Lab/Yun project/duke_result/K562-results-julien/gRNA_prediction/data_new/Hichip/k562_H3K27ac_encode_replicated_peaks_narrowpeaks_ENCFF044JNJ.bed", format = "BED",
#                         extraCols = extraCols_narrowPeak)
hic <- read_tsv("/proj/milovelab/mu/dukeproj/K562-results-julien/gRNA_prediction/K562_combined_r1r2r3.10k.2.peaks.bedpe")
library(GenomicRanges)
library(InteractionSet)
# Converting data.frame to GInteraction
convertToGI <- function(df){
  row.regions <- GRanges(df$chr1, IRanges(df$start1,df$end1))# interaction start
  col.regions <- GRanges(df$chr2, IRanges(df$start2,df$end2))# interaction end
  gi <- GInteractions(row.regions, col.regions)
  mcols(gi) <- df[,7:14] # Interaction frequencies
  return(gi)
}

hic.gi <- convertToGI(hic)
hic.gi$log10fdr <- -log10(hic.gi$fdr)
hic.gi$log10fdr[is.infinite(hic.gi$log10fdr)] <- max(hic.gi$log10fdr[is.finite(hic.gi$log10fdr)])
library(readr)
library(plyranges)
setwd("/pine/scr/w/a/wancen/gRNA_prediction/")
# dat = read_csv("./data_new/wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05.csv")
# colnames(dat)
dat <- data_sub
gr <- makeGRangesFromDataFrame(dat, keep.extra.columns=TRUE,strand.field = "strand.y",
                                start.field="chromStart",end.field="chromEnd")
gr <- makeGRangesFromDataFrame(dat, keep.extra.columns=TRUE,strand.field = "strand",
                                start.field="start",end.field="end")

olap <- findOverlaps(gr,hic.gi) %>% as.data.frame()
olap$log10fdr <- hic.gi$log10fdr[olap$subjectHits]
sumlogfdr <- olap %>% group_by(queryHits) %>% summarise(sumlog10fdr=sum(log10fdr))

countolap <- countOverlaps(gr,hic.gi)
ref = read_tsv("/proj/milovelab/mu/dukeproj/K562-results-julien/gRNA_prediction/gene_info.gencode.v28lift37.txt")
refgr <- makeGRangesFromDataFrame(ref, keep.extra.columns=TRUE,
                               start.field="Promoter_Start",end.field="Promoter_End")

out <- linkOverlaps(hic.gi, gr, refgr) %>% as.data.frame()
out <- unique(out[,c(1,2)])
out$log10fdr <- hic.gi$log10fdr[out$query]
prom <- out %>% dplyr::count(subject1)
promfdr <- out %>% group_by(subject1) %>% summarise(sumlog10fdr=sum(log10fdr))

dat<-data_sub
dat$promnumber <- rep(0,nrow(dat))
dat$promnumber[prom$subject1]<-prom$n
dat$promlog10fdr <- rep(0,nrow(dat))
dat$promlog10fdr[promfdr$subject1]<-promfdr$sumlog10fdr
dat$olapnumber <- countolap
dat$sumlog10fdr <- rep(0,nrow(dat))
dat$sumlog10fdr[sumlogfdr$queryHits]<-sumlogfdr$sumlog10fdr
write.table(dat, file = "/proj/milovelab/mu/dukeproj/K562-supp-tables/screen_hg19_uniqe.txt", sep = "\t",quote=FALSE,
            row.names = FALSE, col.names = TRUE)

## read H3K27ac data and change it to GRanges

H3K27ac_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/H3K27ac_cpm1kb.csv",header = T)
H3K27ac_cpm1kb_gr <- makeGRangesFromDataFrame(H3K27ac_cpm1kb, keep.extra.columns=TRUE)

## read DNAse data
DNAse_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/DNase_cpm1kb.csv",header = T)
DNAse_cpm1kb_gr <- makeGRangesFromDataFrame(DNAse_cpm1kb, keep.extra.columns=TRUE)

## read ATAC-seq data
ATAC_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/ATAC_cpm1kb.csv",header = T)
ATAC_cpm1kb_gr <- makeGRangesFromDataFrame(ATAC_cpm1kb, keep.extra.columns=TRUE)

## read H3Kme4 data
H3Kme4_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/H3Kme4_cpm1kb.csv",header = T)
H3Kme4_cpm1kb_gr <- makeGRangesFromDataFrame(H3Kme4_cpm1kb, keep.extra.columns=TRUE)

## read TF_GATA2 data
TF_GATA2_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/TF_GATA2_cpm1kb.csv",header = T)
TF_GATA2_cpm1kb_gr <- makeGRangesFromDataFrame(TF_GATA2_cpm1kb, keep.extra.columns=TRUE)

## read TF_TAL1 data
TF_TAL1_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/TF_TAL1_cpm1kb.csv",header = T)
TF_TAL1_cpm1kb_gr <- makeGRangesFromDataFrame(TF_TAL1_cpm1kb, keep.extra.columns=TRUE)

## read TF_GATA2 data
TF_MYC_cpm1kb<-read.table("/proj/milovelab/mu/dukeproj/K562-supp-tables/TF_MYC_cpm1kb.csv",header = T)
TF_MYC_cpm1kb_gr <- makeGRangesFromDataFrame(TF_MYC_cpm1kb, keep.extra.columns=TRUE)

library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
H3K27ac_cpm1kb_gr = liftOver(H3K27ac_cpm1kb_gr, ch)
H3K27ac_cpm1kb_gr = unlist(H3K27ac_cpm1kb_gr)
genome(H3K27ac_cpm1kb_gr) = "hg19"

DNAse_cpm1kb_gr = liftOver(DNAse_cpm1kb_gr, ch)
DNAse_cpm1kb_gr = unlist(DNAse_cpm1kb_gr)
genome(DNAse_cpm1kb_gr) = "hg19"

ATAC_cpm1kb_gr = liftOver(ATAC_cpm1kb_gr, ch)
ATAC_cpm1kb_gr = unlist(ATAC_cpm1kb_gr)
genome(ATAC_cpm1kb_gr) = "hg19"

H3Kme4_cpm1kb_gr = liftOver(H3Kme4_cpm1kb_gr, ch)
H3Kme4_cpm1kb_gr = unlist(H3Kme4_cpm1kb_gr)
genome(H3Kme4_cpm1kb_gr) = "hg19"

TF_GATA2_cpm1kb_gr = liftOver(TF_GATA2_cpm1kb_gr, ch)
TF_GATA2_cpm1kb_gr = unlist(TF_GATA2_cpm1kb_gr)
genome(TF_GATA2_cpm1kb_gr) = "hg19"

TF_TAL1_cpm1kb_gr = liftOver(TF_TAL1_cpm1kb_gr, ch)
TF_TAL1_cpm1kb_gr = unlist(TF_TAL1_cpm1kb_gr)
genome(TF_TAL1_cpm1kb_gr) = "hg19"

TF_MYC_cpm1kb_gr = liftOver(TF_MYC_cpm1kb_gr, ch)
TF_MYC_cpm1kb_gr = unlist(TF_MYC_cpm1kb_gr)
genome(TF_MYC_cpm1kb_gr) = "hg19"

# table5 <- read_table2(file="/proj/milovelab/mu/dukeproj/K562-supp-tables/supplementary_table_5_DHS_summary_results.csv")
table5 <- read_table2("/proj/milovelab/mu/dukeproj/K562-supp-tables/supplementary_table_5_DHS_summary_results.txt.gz")
table5_gr <- makeGRangesFromDataFrame(table5, keep.extra.columns=TRUE,
                                      start.field="chromStart",end.field="chromEnd")

olap <- findOverlaps(table5_gr,H3K27ac_cpm1kb_gr) %>% as.data.frame()
olap$H3k27ac_CPM_1Kb_new <- H3K27ac_cpm1kb_gr$cpm1kb[olap$subjectHits]
H3k27ac_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(H3k27ac_CPM_1Kb_new=sum(H3k27ac_CPM_1Kb_new))
table5$H3k27ac_CPM_1Kb_new <-0
table5$H3k27ac_CPM_1Kb_new[H3k27ac_CPM_1Kb_new$queryHits]<-H3k27ac_CPM_1Kb_new$H3k27ac_CPM_1Kb_new

olap <- findOverlaps(table5_gr,DNAse_cpm1kb_gr) %>% as.data.frame()
olap$DNAse_CPM_1Kb_new <- DNAse_cpm1kb_gr$cpm1kb[olap$subjectHits]
DNAse_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(DNAse_CPM_1Kb_new=sum(DNAse_CPM_1Kb_new))
table5$DNAse_CPM_1Kb_new <-0
table5$DNAse_CPM_1Kb_new[DNAse_CPM_1Kb_new$queryHits]<-DNAse_CPM_1Kb_new$DNAse_CPM_1Kb_new

olap <- findOverlaps(table5_gr,ATAC_cpm1kb_gr) %>% as.data.frame()
olap$ATAC_CPM_1Kb_new <- ATAC_cpm1kb_gr$cpm1kb[olap$subjectHits]
ATAC_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(ATAC_CPM_1Kb_new=sum(ATAC_CPM_1Kb_new))
table5$ATAC_CPM_1Kb_new <-0
table5$ATAC_CPM_1Kb_new[ATAC_CPM_1Kb_new$queryHits]<-ATAC_CPM_1Kb_new$ATAC_CPM_1Kb_new

olap <- findOverlaps(table5_gr,H3Kme4_cpm1kb_gr) %>% as.data.frame()
olap$H3Kme4_CPM_1Kb_new <- H3Kme4_cpm1kb_gr$cpm1kb[olap$subjectHits]
H3Kme4_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(H3Kme4_CPM_1Kb_new=sum(H3Kme4_CPM_1Kb_new))
table5$H3Kme4_CPM_1Kb_new <-0
table5$H3Kme4_CPM_1Kb_new[H3Kme4_CPM_1Kb_new$queryHits]<-H3Kme4_CPM_1Kb_new$H3Kme4_CPM_1Kb_new

olap <- findOverlaps(table5_gr,TF_GATA2_cpm1kb_gr) %>% as.data.frame()
olap$TF_GATA2_CPM_1Kb_new <- TF_GATA2_cpm1kb_gr$cpm1kb[olap$subjectHits]
TF_GATA2_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(TF_GATA2_CPM_1Kb_new=sum(TF_GATA2_CPM_1Kb_new))
table5$TF_GATA2_CPM_1Kb_new <-0
table5$TF_GATA2_CPM_1Kb_new[TF_GATA2_CPM_1Kb_new$queryHits]<-TF_GATA2_CPM_1Kb_new$TF_GATA2_CPM_1Kb_new

olap <- findOverlaps(table5_gr,TF_TAL1_cpm1kb_gr) %>% as.data.frame()
olap$TF_TAL1_CPM_1Kb_new <- TF_TAL1_cpm1kb_gr$cpm1kb[olap$subjectHits]
TF_TAL1_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(TF_TAL1_CPM_1Kb_new=sum(TF_TAL1_CPM_1Kb_new))
table5$TF_TAL1_CPM_1Kb_new <-0
table5$TF_TAL1_CPM_1Kb_new[TF_TAL1_CPM_1Kb_new$queryHits]<-TF_TAL1_CPM_1Kb_new$TF_TAL1_CPM_1Kb_new

olap <- findOverlaps(table5_gr,TF_MYC_cpm1kb_gr) %>% as.data.frame()
olap$TF_MYC_CPM_1Kb_new <- TF_MYC_cpm1kb_gr$cpm1kb[olap$subjectHits]
TF_MYC_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(TF_MYC_CPM_1Kb_new=sum(TF_MYC_CPM_1Kb_new))
table5$TF_MYC_CPM_1Kb_new <-0
table5$TF_MYC_CPM_1Kb_new[TF_MYC_CPM_1Kb_new$queryHits]<-TF_MYC_CPM_1Kb_new$TF_MYC_CPM_1Kb_new

# write.table(table5,file="/proj/milovelab/mu/dukeproj/K562-supp-tables/supplementary_table_5_DHS_summary_results.csv",row.names = F)

## add H3K27ac to unique table1
dat<- read_table2("/proj/milovelab/mu/dukeproj/K562-supp-tables/screen_hg19_uniqe.txt")

dat1 = read_csv("/proj/milovelab/mu/dukeproj/K562-supp-tables/discovery_screen_k562_sgrna_deseq2_results_hg19_new.csv", guess_max = 10000)
dat1 = dat1[!duplicated(dat1$protospacer),]
dat<-dat %>% left_join(dat1[,7:ncol(dat1)],by="protospacer")

table1 <- read_table2("/proj/milovelab/mu/dukeproj/K562-results-julien/wgCERES-gRNAs.annotated.txt.gz")

table5_select <- table5[,c("name","annotation_wg","H3k27ac_CPM_1Kb_new", "DNAse_CPM_1Kb_new", "ATAC_CPM_1Kb_new",
                           "H3Kme4_CPM_1Kb_new","TF_GATA2_CPM_1Kb_new", "TF_TAL1_CPM_1Kb_new","TF_MYC_CPM_1Kb_new")]
install.packages("picante")
library(picante)
#Insert the name of your dataset in the code below
cor <- cor.table(table5_select[,3:9], cor.method="pearson")

table5_select$annotation_wg_orig = table5_select$annotation_wg
table5_select$annotation_wg <- ifelse(table5_select$annotation_wg == "Promoter", "Promoter","enhancer")

table1_select <- table1[,-c(1:6,8)]

data <- dat %>% left_join(table1_select,by="gRNAid")
data <- data %>% left_join(table5_select, by=c("DHS"="name"))
hist(log(data$H3k27ac_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(H3K27ac_cpm_per1kb)",main = "promoter")
hist(log(data$H3k27ac_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),xlab = "log(H3K27ac_cpm_per1kb)",main = "enhancer")
hist(log(data$H3k27ac_CPM_1Kb_new[which(data$annotation_wg_orig=="Intron")]+1),xlab = "log(H3K27ac_cpm_per1kb)",main = "enhancer")

hist(log(data$DNase_CPM_per_1kbp[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(DNAse_cpm_per1kb)",main = "promoter")
hist(log(data$DNase_CPM_per_1kbp[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(DNAse_cpm_per1kb)",main = "enhancer")

hist(log(data$ATAC_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(ATAC_cpm_per1kb)",main = "promoter")
hist(log(data$ATAC_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(ATAC_cpm_per1kb)",main = "enhancer")

hist(log(data$H3Kme4_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(H3Kme4_cpm_per1kb)",main = "promoter")
hist(log(data$H3Kme4_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(H3Kme4_cpm_per1kb)",main = "enhancer")

hist(log(data$TF_GATA2_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(TF_GATA2_cpm_per1kb)",main = "promoter")
hist(log(data$TF_GATA2_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(TF_GATA2_cpm_per1kb)",main = "enhancer")

hist(log(data$TF_TAL1_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(TF_TAL1_cpm_per1kb)",main = "promoter")
hist(log(data$TF_TAL1_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(TF_TAL1_cpm_per1kb)",main = "enhancer")

hist(log(data$TF_MYC_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),breaks=40,xlab = "log(TF_MYC_cpm_per1kb)",main = "promoter")
hist(log(data$TF_MYC_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),breaks=40,xlab = "log(TF_MYC_cpm_per1kb)",main = "enhancer")

plot(x=log(data$H3k27ac_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="Promoter")],
     xlab = "log(H3K27ac_cpm_per1kb)",main = "promoter",ylab="logFC")
plot(x=log(data$H3k27ac_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="enhancer")],
     xlab = "log(H3K27ac_cpm_per1kb)",main = "enhancer",ylab="logFC")

plot(x=log(data$DNAse_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="Promoter")],
     xlab = "log(DNAse_cpm_per1kb)",main = "promoter",ylab="logFC")
plot(x=log(data$DNAse_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="enhancer")],
     xlab = "log(DNAse_cpm_per1kb)",main = "enhancer",ylab="logFC")

plot(x=log(data$ATAC_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="Promoter")],
     xlab = "log(ATAC_cpm_per1kb)",main = "promoter",ylab="logFC")
plot(x=log(data$ATAC_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="enhancer")],
     xlab = "log(ATAC_cpm_per1kb)",main = "enhancer",ylab="logFC")

plot(x=log(data$TF_GATA2_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="Promoter")],
     xlab = "log(TF_GATA2_cpm_per1kb)",main = "promoter",ylab="logFC")
plot(x=log(data$TF_TAL1_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="enhancer")],
     xlab = "log(TF_TAL1_cpm_per1kb)",main = "enhancer",ylab="logFC")

plot(x=log(data$TF_MYC_CPM_1Kb_new[which(data$annotation_wg=="Promoter")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="Promoter")],
     xlab = "log(TF_GATA2_cpm_per1kb)",main = "promoter",ylab="logFC")
plot(x=log(data$TF_MYC_CPM_1Kb_new[which(data$annotation_wg=="enhancer")]+1),
     y=data$log2FoldChange[which(data$annotation_wg=="enhancer")],
     xlab = "log(TF_TAL1_cpm_per1kb)",main = "enhancer",ylab="logFC")
# protospacer <- dat[,"protospacer"]
# a <- unique(sapply(protospacer, nchar))
# position <- matrix(NA,ncol=a,nrow=nrow(dat)) %>% as.data.frame
# for(i in 1:a){
#   position[,i] <- as.factor(substr(protospacer$protospacer,i,i))
# }
# colnames(position) <- paste("position",seq_len(20),sep = "_")
# position2 <- matrix(NA,ncol=19,nrow=nrow(dat)) %>% as.data.frame
# for(i in 1:19){
#   position2[,i] <- as.factor(substr(protospacer$protospacer,i,i+1))
# }
# colnames(position2) <- paste("position",seq_len(19),sep = "_")
# 
# library(caret)
# f <- paste("~", paste(colnames(position), collapse="+"))
# dummy <- dummyVars(f, data=position,sep = "_")
# newdata <- data.frame(predict(dummy, newdata = position))
# 
# f <- paste("~", paste(colnames(position2), collapse="+"))
# dummy <- dummyVars(f, data=position2,sep = "_")
# newdata2 <- data.frame(predict(dummy, newdata = position2))
# 
# dinucleotide<- levels(position2$position_1)
# count <- matrix(NA,ncol=16,nrow=nrow(dat)) %>% as.data.frame
# for(i in 1:16){
#   count[,i] <- str_count(protospacer$protospacer,dinucleotide[i])
# }
# colnames(count) <- paste0(dinucleotide,"count")

## select promoter
data$vc_sqrt_sum <- as.numeric(as.character(data$vc_sqrt_sum))
dat_pro = data %>% filter(annotation_wg=="Promoter")
dat_sel = dat_pro %>%
  filter(pvalue <= 0.05)
dat_sel = dat_sel %>%
  mutate(significant = ifelse(padj <= 0.05, 1, 0))
dat_nonsig = dat_pro %>%
  filter(pvalue > 0.05)

write_csv(dat_sel, "./data_new/wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05.csv")
dat_sel = read_csv("./data_new/wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05.csv")
library(tidyverse)
dat_sel<-dat_sel %>% left_join(data_range[,6:ncol(data_range)],by=c("gRNAid"="GeneID"))
dat_sel = dat_sel %>%
  mutate(direction_sig = ifelse(padj > 0.05, 0, ifelse(log2FoldChange<=0,1,2)))
dat_sel$vc_sqrt_sum <- as.numeric(as.character(dat_sel$vc_sqrt_sum))
write_csv(dat_nonsig, "./data_new/wgCERES-gRNAs-k562-discovery-screen-pro_rawp-nonsignificant.csv")

library(caret)
set.seed(2021)
dat_padj_sel = dat_sel %>%
  filter(padj <= 0.2)

trainid = createDataPartition(dat_padj_sel$significant, p=0.8, list=FALSE)
trainid = createDataPartition(dat_padj_sel$direction_sig, p=0.8, list=FALSE)
train = dat_padj_sel[trainid[,1],]
test = dat_padj_sel[-trainid[,1],]

write_csv(train, "./data_new/wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05-padj0.2-train2.csv")
write_csv(test, "./data_new/wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05-padj0.2-test2.csv")


##### binary outcome prediction: filter out padj between 0.05 and 0.2
dat_sel_binary = dat_sel %>%
  filter(padj<=0.05 | padj >= 0.2)

set.seed(2021)
trainid <- createDataPartition(dat_sel_binary$significant, p=0.8, list=FALSE)
trainid = createDataPartition(dat_sel_binary$direction_sig, p=0.8, list=FALSE)
train = dat_sel_binary[trainid[,1],]
test = dat_sel_binary[-trainid[,1],]

write_csv(train, "./data_new/wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05-binary-train2.csv")
write_csv(test, "./data_new/wgCERES-gRNAs-k562-discovery-screen-pro_rawp0.05-binary-test2.csv")

## select enhancer
dat_enh = data %>% filter(annotation_wg=="enhancer")
dat_sel = dat_enh %>%
  filter(pvalue <= 0.05)
dat_sel = dat_sel %>%
  mutate(significant = ifelse(padj <= 0.05, 1, 0))
dat_sel$vc_sqrt_sum <- as.numeric(as.character(dat_sel$vc_sqrt_sum))
dat_sel = dat_sel %>%
  mutate(direction_sig = ifelse(padj > 0.05, 0, ifelse(log2FoldChange<=0,1,2)))
dat_nonsig = dat_enh %>%
  filter(pvalue > 0.05)

write_csv(dat_sel, "./data_new/wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05.csv")
write_csv(dat_nonsig, "./data_new/wgCERES-gRNAs-k562-discovery-screen-enh_rawp-nonsignificant.csv")

library(caret)
set.seed(2021)
dat_padj_sel = dat_sel %>%
  filter(padj <= 0.2)

trainid = createDataPartition(dat_padj_sel$significant, p=0.8, list=FALSE)
trainid = createDataPartition(dat_padj_sel$direction_sig, p=0.8, list=FALSE)
train = dat_padj_sel[trainid[,1],]
test = dat_padj_sel[-trainid[,1],]

write_csv(train, "./data_new/wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-padj0.2-train.csv")
write_csv(test, "./data_new/wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-padj0.2-test.csv")


##### binary outcome prediction: filter out padj between 0.05 and 0.2
dat_sel_binary = dat_sel %>%
  filter(padj<=0.05 | padj >= 0.2)

set.seed(2021)
trainid <- createDataPartition(dat_sel_binary$significant, p=0.8, list=FALSE)
trainid = createDataPartition(dat_sel_binary$direction_sig, p=0.8, list=FALSE)
train = dat_sel_binary[trainid[,1],]
test = dat_sel_binary[-trainid[,1],]

write_csv(train, "./data_new/wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-binary-train.csv")
write_csv(test, "./data_new/wgCERES-gRNAs-k562-discovery-screen-enh_rawp0.05-binary-test.csv")

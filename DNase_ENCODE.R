## two replicates
file1 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF425WDA.bam"
file2 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF205FNC.bam"
library(Rsamtools)
indexBam(file1)
indexBam(file2)
library(GenomicRanges)
library(GenomeInfoDb)
si <- Seqinfo(genome="hg38")
si <- keepSeqlevels(si, value=paste0("chr",c(1:22,"X")))
tiles <- tileGenome(si, tilewidth=1e3, cut.last.tile.in.chrom=TRUE)
tiles_df <- as.data.frame(tiles)
tiles_df <- data.frame(GeneID = 1:nrow(tiles_df),
                       tiles_df[,c(1,2,3,5)])
names(tiles_df)[2:5] <- c("Chr","Start","End","Strand")

library(Rsubread)
res1 <- featureCounts(files=file1,annot.inbuilt = "hg38", annot.ext=data_range)
res2 <- featureCounts(files=file2,annot.inbuilt = "hg38", annot.ext=data_range,isPairedEnd = TRUE)
hist(log10(res1$counts+1))
hist(log10(res2$counts+1))
avg_count <- (res1$counts+res2$counts)/2
data_range$DNase_cpm_slide1000 <- avg_count*1e6/sum(avg_count)

write.table(tiles_df,file="/proj/milovelab/mu/dukeproj/K562-supp-tables/DNase_cpm1kb.csv",row.names = F)
write.table(data_range,file="/proj/milovelab/mu/dukeproj/K562-supp-tables/DNase_cpm1kb_slide.csv",row.names = F)
## ATACseq #############################################
########################################################
## two replicates
file1 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF512VEZ.bam"
file2 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF987XOV.bam"
library(Rsamtools)
indexBam(file1)
indexBam(file2)
library(GenomicRanges)
library(GenomeInfoDb)
si <- Seqinfo(genome="hg38")
si <- keepSeqlevels(si, value=paste0("chr",c(1:22,"X")))
tiles <- tileGenome(si, tilewidth=1e3, cut.last.tile.in.chrom=TRUE)
tiles_df <- as.data.frame(tiles)
tiles_df <- data.frame(GeneID = 1:nrow(tiles_df),
                       tiles_df[,c(1,2,3,5)])
names(tiles_df)[2:5] <- c("Chr","Start","End","Strand")

library(Rsubread)
res1 <- featureCounts(files=file1,annot.inbuilt = "hg38", annot.ext=tiles_df,isPairedEnd = TRUE)
res2 <- featureCounts(files=file2,annot.inbuilt = "hg38", annot.ext=tiles_df,isPairedEnd = TRUE)
hist(log10(res1$counts+1))
hist(log10(res2$counts+1))
avg_count <- (res1$counts+res2$counts)/2
tiles_df$cpm1kb <- avg_count*1e6/sum(avg_count)

write.table(tiles_df,file="/proj/milovelab/mu/dukeproj/K562-supp-tables/ATAC_cpm1kb.csv",row.names = F)


## H3Kme4 #############################################
########################################################
## two replicates
file1 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF656DMV.bam"
file2 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF440ARP.bam"
library(Rsamtools)
indexBam(file1)
indexBam(file2)
library(GenomicRanges)
library(GenomeInfoDb)
si <- Seqinfo(genome="hg38")
si <- keepSeqlevels(si, value=paste0("chr",c(1:22,"X")))
tiles <- tileGenome(si, tilewidth=1e3, cut.last.tile.in.chrom=TRUE)
tiles_df <- as.data.frame(tiles)
tiles_df <- data.frame(GeneID = 1:nrow(tiles_df),
                       tiles_df[,c(1,2,3,5)])
names(tiles_df)[2:5] <- c("Chr","Start","End","Strand")

library(Rsubread)
res1 <- featureCounts(files=file1,annot.inbuilt = "hg38", annot.ext=tiles_df)
res2 <- featureCounts(files=file2,annot.inbuilt = "hg38", annot.ext=tiles_df)
hist(log10(res1$counts+1))
hist(log10(res2$counts+1))
avg_count <- (res1$counts %>% unname()+res2$counts %>% unname())/2
tiles_df$cpm1kb <- avg_count*1e6/sum(avg_count)

write.table(tiles_df,file="/proj/milovelab/mu/dukeproj/K562-supp-tables/H3Kme4_cpm1kb.csv",row.names = F)

## TF Chip-seq(GATA2) #############################################
########################################################
## two replicates
file1 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF307XGY.bam"
file2 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF032UMK.bam"
library(Rsamtools)
indexBam(file1)
indexBam(file2)
library(GenomicRanges)
library(GenomeInfoDb)
si <- Seqinfo(genome="hg38")
si <- keepSeqlevels(si, value=paste0("chr",c(1:22,"X")))
tiles <- tileGenome(si, tilewidth=1e3, cut.last.tile.in.chrom=TRUE)
tiles_df <- as.data.frame(tiles)
tiles_df <- data.frame(GeneID = 1:nrow(tiles_df),
                       tiles_df[,c(1,2,3,5)])
names(tiles_df)[2:5] <- c("Chr","Start","End","Strand")

library(Rsubread)
res1 <- featureCounts(files=file1,annot.inbuilt = "hg38", annot.ext=tiles_df)
res2 <- featureCounts(files=file2,annot.inbuilt = "hg38", annot.ext=tiles_df)
hist(log10(res1$counts+1))
hist(log10(res2$counts+1))
avg_count <- (res1$counts %>% unname()+res2$counts %>% unname())/2
tiles_df$cpm1kb <- avg_count*1e6/sum(avg_count)

write.table(tiles_df,file="/proj/milovelab/mu/dukeproj/K562-supp-tables/TF_GATA2_cpm1kb.csv",row.names = F)

## TF Chip-seq(TAL1) #############################################
########################################################
## two replicates
file1 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF304UTB.bam"
file2 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF631RZK.bam"
library(Rsamtools)
indexBam(file1)
indexBam(file2)
library(GenomicRanges)
library(GenomeInfoDb)
si <- Seqinfo(genome="hg38")
si <- keepSeqlevels(si, value=paste0("chr",c(1:22,"X")))
tiles <- tileGenome(si, tilewidth=1e3, cut.last.tile.in.chrom=TRUE)
tiles_df <- as.data.frame(tiles)
tiles_df <- data.frame(GeneID = 1:nrow(tiles_df),
                       tiles_df[,c(1,2,3,5)])
names(tiles_df)[2:5] <- c("Chr","Start","End","Strand")

library(Rsubread)
res1 <- featureCounts(files=file1,annot.inbuilt = "hg38", annot.ext=tiles_df,isPairedEnd = TRUE)
res2 <- featureCounts(files=file2,annot.inbuilt = "hg38", annot.ext=tiles_df,isPairedEnd = TRUE)
hist(log10(res1$counts+1))
hist(log10(res2$counts+1))
avg_count <- (res1$counts %>% unname()+res2$counts %>% unname())/2
tiles_df$cpm1kb <- avg_count*1e6/sum(avg_count)

write.table(tiles_df,file="/proj/milovelab/mu/dukeproj/K562-supp-tables/TF_TAL1_cpm1kb.csv",row.names = F)

## TF Chip-seq(MYC) #############################################
########################################################
## two replicates
file1 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF503ZCR.bam"
file2 <- "/pine/scr/w/a/wancen/gRNA_prediction/H3K27/ENCFF604CDP.bam"
library(Rsamtools)
indexBam(file1)
indexBam(file2)
library(GenomicRanges)
library(GenomeInfoDb)
si <- Seqinfo(genome="hg38")
si <- keepSeqlevels(si, value=paste0("chr",c(1:22,"X")))
tiles <- tileGenome(si, tilewidth=1e3, cut.last.tile.in.chrom=TRUE)
tiles_df <- as.data.frame(tiles)
tiles_df <- data.frame(GeneID = 1:nrow(tiles_df),
                       tiles_df[,c(1,2,3,5)])
names(tiles_df)[2:5] <- c("Chr","Start","End","Strand")

library(Rsubread)
res1 <- featureCounts(files=file1,annot.inbuilt = "hg38", annot.ext=tiles_df,isPairedEnd = TRUE)
res2 <- featureCounts(files=file2,annot.inbuilt = "hg38", annot.ext=tiles_df,isPairedEnd = TRUE)
hist(log10(res1$counts+1))
hist(log10(res2$counts+1))
avg_count <- (res1$counts %>% unname()+res2$counts %>% unname())/2
tiles_df$cpm1kb <- avg_count*1e6/sum(avg_count)

write.table(tiles_df,file="/proj/milovelab/mu/dukeproj/K562-supp-tables/TF_MYC_cpm1kb.csv",row.names = F)

library(rtracklayer)
library(GenomicRanges)
library(InteractionSet)
library(tidyverse)

hic_pp <- read_tsv("/proj/yunligrp/users/tianyou/gRNA/NPC/source_file/pcHiC/NPC.pp.txt")
hic_po <- read_tsv("/proj/yunligrp/users/tianyou/gRNA/NPC/source_file/pcHiC/NPC.po.txt")
hic = bind_rows(hic_pp, hic_po) %>%
  filter(`-log10(result)` >= 2)

## prepare the coordinate file for liftover to hg38
# hic %>%
#   select(frag1, frag2) %>%
#   mutate(ID = row_number()) %>%
#   pivot_longer(cols = c("frag1", "frag2"), values_to = "frag") %>%
#   separate(frag, into = c("chrom", "bp"), sep = ":") %>%
#   separate(bp, into = c("chromStart", "chromEnd"), sep = "-") %>%
#   mutate(name = paste0(ID,"_",name)) %>%
#   mutate(chromStart = as.numeric(chromStart) - 1) %>%
#   select(chrom, chromStart, chromEnd, name) %>%
#   write_tsv("HiC_position_hg19.bed", col_names = FALSE)

## read in liftovered hg38 locations and original locations
coord_hg19 = read_table("/proj/yunligrp/users/tianyou/gRNA/NPC/HiC_position_hg19.bed", col_names = F)
coord_hg38 = read_table("/proj/yunligrp/users/tianyou/gRNA/NPC/HiC_position_lifthg38.bed", col_names = F)
colnames(coord_hg19) = c("chrom", "chromStart", "chromEnd", "name")
colnames(coord_hg38) = c("chrom", "chromStart", "chromEnd", "name", "score")
coord_hg19 = coord_hg19 %>% 
  mutate(frag = str_c(chrom,str_c(chromStart+1,chromEnd,sep="-"), sep = ":")) %>%
  select(chrom, name, frag)
coord_hg38 = coord_hg38 %>% 
  mutate(frag_hg38 = str_c(chrom,str_c(chromStart+1,chromEnd,sep="-"), sep = ":")) %>%
  select(chrom, name, frag_hg38)

coord_update_lst = left_join(coord_hg19, coord_hg38, by = c("chrom", "name"))
coord_dup = coord_update_lst %>% group_by(name) %>% summarise(count = n()) %>% filter(count > 1)
coord_update_lst = coord_update_lst %>%
  filter(!(name %in% coord_dup$name)) %>%
  separate(name, into = c("rownum", "name"), sep = "_") %>%
  pivot_wider(names_from = name, values_from = c("frag","frag_hg38")) %>%
  select(-rownum) %>%
  rename("frag1" = "frag_frag1", "frag2" = "frag_frag2", "frag1_hg38" = "frag_hg38_frag1",
         "frag2_hg38" = "frag_hg38_frag2")

hic_hg38 = hic %>% inner_join(coord_update_lst, by = c("frag1","frag2")) %>%
  filter(!is.na(frag1_hg38), !is.na(frag2_hg38))

## Separate the frag columns to get start/end position
hic_hg38 = hic_hg38 %>%
  separate(frag1_hg38, into = c("chr1", "bp1"), sep = ":") %>%
  separate(bp1, into = c("start1", "end1"), sep = "-") %>%
  separate(frag2_hg38, into = c("chr2", "bp2"), sep = ":") %>%
  separate(bp2, into = c("start2", "end2"), sep = "-")

hic_hg38$start1 = as.numeric(hic_hg38$start1)
hic_hg38$end1 = as.numeric(hic_hg38$end1)
hic_hg38$start2 = as.numeric(hic_hg38$start2)
hic_hg38$end2 = as.numeric(hic_hg38$end2)


# Converting data.frame to GInteraction
convertToGI <- function(df){
  row.regions <- GRanges(df$chr1, IRanges(df$start1,df$end1))# interaction start
  col.regions <- GRanges(df$chr2, IRanges(df$start2,df$end2))# interaction end
  gi <- GInteractions(row.regions, col.regions)
  mcols(gi) <- df[,10] # Interaction frequencies ## metadata column
  return(gi)
}

hic.gi <- convertToGI(hic_hg38)
library(vroom)
library(plyranges)
library(data.table)
dat <- vroom("/proj/yunligrp/users/tianyou/gRNA/from_Alex/MHC_CDR_adjusted/npc.latentvar.MAST_expect_cells.grna_de.markers.aggr.txt.gz")
dat2 <- vroom("/proj/yunligrp/users/tianyou/gRNA/from_Alex/archive/MHC/npc.expect_cells.grna.de.markers.MAST.annotatedfull.final.update20220117.LRB.tsv")

dat3 <- dat %>% group_by(grna) %>% slice_min(p_val,n=1,with_ties = FALSE) %>% 
  left_join(dat2[,c(7:20,40:47)], by = c("grna" = "protospacer", "gene_symbol"))
dat3 <- dat3[!duplicated(dat3), ] %>% filter(type == "targeting")  ## what caused the duplication?

gr <- makeGRangesFromDataFrame(dat3, keep.extra.columns=TRUE,
                               seqnames.field="grna.chr",
                               strand.field = "grna.strand",
                               start.field="grna.start",end.field="grna.end")

olap <- findOverlaps(gr,hic.gi) %>% as.data.frame()
# olap$log10fdr <- hic.gi$log10fdr[olap$subjectHits]
# sumlogfdr <- olap %>% group_by(queryHits) %>% summarise(sumlog10fdr=sum(log10fdr))

countolap <- countOverlaps(gr,hic.gi)
ref = read_tsv("/proj/milovelab/mu/dukeproj/K562-results-julien/gRNA_prediction/gene_info.gencode.v28lift37.txt")
refgr <- makeGRangesFromDataFrame(ref, keep.extra.columns=TRUE,
                                  start.field="Promoter_Start",end.field="Promoter_End")

out <- linkOverlaps(hic.gi, gr, refgr) %>% as.data.frame()
out <- unique(out[,c(1,2)])
# out$log10fdr <- hic.gi$log10fdr[out$query]
prom <- out %>% dplyr::count(subject1)
# promfdr <- out %>% group_by(subject1) %>% summarise(sumlog10fdr=sum(log10fdr))

dat<-dat3
dat$promnumber <- rep(0,nrow(dat))
dat$promnumber[prom$subject1]<-prom$n
# dat$promlog10fdr <- rep(0,nrow(dat))
# dat$promlog10fdr[promfdr$subject1]<-promfdr$sumlog10fdr
dat$olapnumber <- countolap
# dat$sumlog10fdr <- rep(0,nrow(dat))
# dat$sumlog10fdr[sumlogfdr$queryHits]<-sumlogfdr$sumlog10fdr

## 581 dhs
unique(dat$dhs) %>% length()
dhs_gr <- makeGRangesFromDataFrame(dat, keep.extra.columns=TRUE,
                                   seqnames.field="dhs.chr",
                                   strand.field = "grna.strand",
                                   start.field="dhs.start",end.field="dhs.end")
## read H3K27ac data and change it to GRanges

H3K27ac_cpm1kb<-read.table("/proj/yunligrp/users/tianyou/gRNA/NPC/NPC_H3K27ac_cpm1kb.csv",header = T)
H3K27ac_cpm1kb_gr <- makeGRangesFromDataFrame(H3K27ac_cpm1kb, keep.extra.columns=TRUE)

## read ATAC-seq data
ATAC_cpm1kb<-read.table("/proj/yunligrp/users/tianyou/gRNA/NPC/NPC_ATAC_cpm1kb.csv",header = T)
ATAC_cpm1kb_gr <- makeGRangesFromDataFrame(ATAC_cpm1kb, keep.extra.columns=TRUE)

## read H3K4me3 data
H3K4me3_cpm1kb<-read.table("/proj/yunligrp/users/tianyou/gRNA/NPC/NPC_H3K4me3_cpm1kb.csv",header = T)
H3K4me3_cpm1kb_gr <- makeGRangesFromDataFrame(H3K4me3_cpm1kb, keep.extra.columns=TRUE)

## overlap with dhs level
olap <- findOverlaps(dhs_gr,H3K27ac_cpm1kb_gr) %>% as.data.frame()
olap$H3k27ac_CPM_1Kb_new <- H3K27ac_cpm1kb_gr$cpm1kb[olap$subjectHits]
H3k27ac_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(H3k27ac_CPM_1Kb_new=mean(H3k27ac_CPM_1Kb_new))
dat$H3k27ac_CPM_1Kb_new <-0
dat$H3k27ac_CPM_1Kb_new[H3k27ac_CPM_1Kb_new$queryHits]<-H3k27ac_CPM_1Kb_new$H3k27ac_CPM_1Kb_new

olap <- findOverlaps(dhs_gr,ATAC_cpm1kb_gr) %>% as.data.frame()
olap$ATAC_CPM_1Kb_new <- ATAC_cpm1kb_gr$cpm1kb[olap$subjectHits]
ATAC_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(ATAC_CPM_1Kb_new=mean(ATAC_CPM_1Kb_new))
dat$ATAC_CPM_1Kb_new <-0
dat$ATAC_CPM_1Kb_new[ATAC_CPM_1Kb_new$queryHits]<-ATAC_CPM_1Kb_new$ATAC_CPM_1Kb_new

olap <- findOverlaps(dhs_gr,H3K4me3_cpm1kb_gr) %>% as.data.frame()
olap$H3K4me3_CPM_1Kb_new <- H3K4me3_cpm1kb_gr$cpm1kb[olap$subjectHits]
H3K4me3_CPM_1Kb_new <- olap %>% group_by(queryHits) %>% summarise(H3K4me3_CPM_1Kb_new=mean(H3K4me3_CPM_1Kb_new))
dat$H3K4me3_CPM_1Kb_new <-0
dat$H3K4me3_CPM_1Kb_new[H3K4me3_CPM_1Kb_new$queryHits]<-H3K4me3_CPM_1Kb_new$H3K4me3_CPM_1Kb_new

## deltagb and deltagh
protospacer <- dat$grna
strand = dat$grna.strand
position2 <- matrix(NA,ncol=19,nrow=length(protospacer)) %>% as.data.frame
for(i in 1:19){
  position2[,i] <- as.factor(substr(protospacer,i,i+1))
}
colnames(position2) <- paste("position",seq_len(19),sep = "_")

delta <- c(-1,-2.1,-1.8,-0.9,-0.9,-2.1,-1.7,-0.9,-1.3,-2.7,-2.9,-1.1,-0.6,-1.5,-1.6,-0.2)
names(delta)<-c("TT","GT","CT","AT","TG","GG","CG","AG","TC","GC","CC","AC","TA","GA","CA","AA")
weight= c(1.80, 1.96, 1.90, 2.13, 1.38, 1.46, 1.00, 1.39, 1.51, 1.98, 1.88, 1.72, 2.02, 1.93, 2.08, 1.94, 2.15, 2.04, 2.25)
deltagh<- sapply(1:length(protospacer), function(i){sum(weight * delta[position2[i,] %>% unname() %>% unlist()] )})
dat$deltagh <-deltagh

rna_rna <- c(-0.93,-2.24,-2.08,-1.1,-2.11,-3.26,-2.36,-2.08,-2.35,-3.42,-3.26,-2.24,-1.33,-2.35,-2.11,-0.93)
dna_dna <- c(-1,-1.44,-1.28,-0.88,-1.45,-1.84,-2.17,-1.28,-1.3,-2.24,-1.84,-1.44,-0.58,-1.3,-1.45,-1)
protospacer2 <- ifelse(strand=="+",protospacer,chartr('ATGC', 'TACG', protospacer))
protospacer2 <- ifelse(strand=="+",protospacer2,reverse(protospacer2))

position3 <- matrix(NA,ncol=19,nrow=length(protospacer2)) %>% as.data.frame
for(i in 1:19){
  position3[,i] <- as.factor(substr(protospacer2,i,i+1))
}
colnames(position3) <- paste("position",seq_len(19),sep = "_")
deltagu<- sapply(1:length(protospacer2), function(i){sum(rna_rna[position3[i,] %>% unname() %>% unlist()] )})

deltago<- sapply(1:length(protospacer), function(i){sum(dna_dna[position2[i,] %>% unname() %>% unlist()] )})
deltagb<- deltagh-deltago-deltagu
dat$deltagb <-deltagb
write_csv(dat,file="/proj/yunligrp/users/tianyou/gRNA/NPC/NPC_onegene_full.csv")

dat %>%
  select(grna.chr, grna.start, grna.end, grna, dhs, grna.strand) %>%
  write_tsv("/proj/yunligrp/users/tianyou/gRNA/NPC/NPC_coordinates.tsv", col_names = F)



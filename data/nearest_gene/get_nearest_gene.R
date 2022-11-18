library(tidyverse)
library(data.table)

## Part to obtain the Ensemble v105 data frame
# library(AnnotationHub)
# ah <- AnnotationHub()
# EnsDbv105 = query(ah, c("EnsDb", "v105", "Homo sapiens"))
# EnsDbv105 = EnsDbv105[["AH98047"]]
# tx <- transcripts(EnsDbv105, columns=c('tx_id', 'tx_id_version', 'tx_biotype', 'gene_id', 'gene_name', 'gene_seq_start', 'gene_seq_end'))
# tx = tx %>% as_tibble() %>% dplyr::filter(tx_biotype == "protein_coding")
# tx = tx %>% 
#   dplyr::select(seqnames, start, end, gene_id, gene_name, strand) %>% 
#   dplyr::filter(seqnames %in% as.character(seq(1,22))) %>%
#   mutate(seqnames = str_c("chr", seqnames)) %>%
#   write_tsv("/proj/yunligrp/users/tianyou/gRNA/OGEE/Ensemble105.bed", col_names = FALSE)


dir = "/proj/yunligrp/users/tianyou/gRNA/OGEE"
data.tbl <- fread(file.path(dir, "NPC_hg38_annotated.tsv"), 
                  col.names = c("chr","start","end","grna","dhs","strand",
                                "gchr","gstart","gend","gene","score","gstrand","distance"))

OGEE = fread(file.path(dir, "OGEE_essentiality.txt"))[,-1]
OGEE_merged = data.tbl %>% 
  left_join(OGEE, by = c("gene" = "locus")) %>% 
  mutate(distToTSS = pmin(abs(start - gstart), abs(end - gstart))) %>%
  replace_na(list(OGEE_prop_Essential = 0)) %>%
  select(chr, start, end, grna, dhs, gene, distance, OGEE_prop_Essential, geneId, distToTSS) %>% 
  unique() %>%
  group_by(chr, start, end, grna) %>%
  slice_max(OGEE_prop_Essential) %>% 
  slice_min(distToTSS)

write_tsv(OGEE_merged, file.path(dir, "NPC_OGEE.txt"))


dir = "/proj/yunligrp/users/tianyou/gRNA/OGEE"
data.tbl <- fread(file.path(dir, "NPC_hg38_annotated_downstream.tsv"), 
                  col.names = c("chr","start","end","grna","dhs","strand",
                                "gchr","gstart","gend","gene","score","gstrand","distance"))

OGEE = fread(file.path(dir, "OGEE_essentiality.txt"))[,-1]
OGEE_merged = data.tbl %>% 
  left_join(OGEE, by = c("gene" = "locus")) %>% 
  mutate(distToTSS = pmin(abs(start - gstart), abs(end - gstart))) %>%
  replace_na(list(OGEE_prop_Essential = 0)) %>%
  select(chr, start, end, grna, dhs, gene, distance, OGEE_prop_Essential, geneId, distToTSS) %>% 
  unique() %>%
  group_by(chr, start, end, grna) %>%
  slice_max(OGEE_prop_Essential) %>% 
  slice_min(distToTSS)

write_tsv(OGEE_merged, file.path(dir, "NPC_OGEE_downstream.txt"))

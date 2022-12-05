library(tidyverse)
library(vroom)
library(data.table)

scdir = "/proj/yunligrp/users/tianyou/gRNA/data/single-cell"

dat <- vroom("/proj/milovelab/mu/dukeproj/data/scOct4th/grna.de.markers.nFeature_RNA.with_grna_coords.txt.gz")
dat2 <- vroom("/pine/scr/t/i/tianyou/Patrick/single-cell/single-cell-data.csv")
colnames(dat)
unique(dat$grna) %>% length()  

dat2 <- dat2[,c(7,14:17)] %>% as.data.table() %>% unique()
dat3 <- dat %>% group_by(grna) %>% slice_min(p_val,n=1,with_ties = FALSE) %>% left_join(dat2, by = c("grna"))
dat3 <- dat3 %>% filter(!is.na(dhs_start))

write_tsv(dat3, file.path(scdir, "grna.de.markers.nFeature_RNA.with_grna_coords_dhs.txt"))
dat3 %>%
  select(chrom, start, end, grna, dhs_id, strand) %>%
  write_tsv("/proj/yunligrp/users/tianyou/gRNA/OGEE/sc_coordinates.tsv", col_names = F)


### get nearest gene
dir = "/proj/yunligrp/users/tianyou/gRNA/OGEE"
data.tbl <- fread(file.path(dir, "sc_annotated.tsv"), 
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

write_tsv(OGEE_merged, file.path(dir, "sc_OGEE.txt"))



# Obtaining the information of nearest gene

We obtain the Ensemble 105 dataframe from AnnotationHub (see get_nearest_gene.R).

First sort the reference file and coordinate files:
dir="/proj/yunligrp/users/tianyou/gRNA/OGEE/"
<!-- awk 'BEGIN{OFS="\t"} NR>1 {print $3,$5,$6,$2,$12,$4}' ${dir}/ALLgencodev39_hg38.tsv | sort -k1,1 -k2,2n > ${dir}/ALLgencodev39_hg38_sorted.bed -->
sort -k1,1 -k2,2n ${dir}/Ensemble105.bed > ${dir}/Ensemble105_sorted.bed
sort -k1,1 -k2,2n /proj/yunligrp/users/tianyou/gRNA/NPC/NPC_coordinates.tsv > ${dir}/NPC_hg38_sorted.bed

(1) Does not requiring downstream
bedtools closest -s -D a -a ${dir}/NPC_hg38_sorted.bed -b ${dir}/Ensemble105_sorted.bed > NPC_hg38_annotated.tsv
(2) Requiring downstream
bedtools closest -s -D a -iu -a ${dir}/NPC_hg38_sorted.bed -b ${dir}/Ensemble105_sorted.bed > NPC_hg38_annotated_downstream.tsv

Then use get_nearest_gene.R to merge these two results with OGEE result to create the variable.
# Processing raw ATAC-seq files for NPC
From the following paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6850896/

File location: 

NarrowPeak format: 	https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115046

SRA raw files: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA473806&f=assay_type_s%3An%3Aatac-seq%3Ac&o=acc_s%3Aa

## Step 1: download the data
Select all and download the accession list. Then load the module sratoolkit on longleaf 
and follow the instructions at https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/ to configure SRA toolkit and download public data.
Use script download.sh and also click on the SRA run page to download metadata.
Files are downloaded into: /pine/scr/t/i/tianyou/sra/sra
We only used files corresponding to 7TP 72hr, which is SRR7230429 and SRR7230436.

## Step 2: install the pipeline environment
(1) Follow the instruction file to build the conda environment. The current anaconda version on longleaf could successfully build the environment.

(2) Download the reference genome (hg38)

`bash scripts/download_genome_data.sh hg38 /pine/scr/t/i/tianyou/ATAC-seq/reference_genome_hg38`

(3) Create the json file

By examining the file names of the NarrowPeak files, we could see which experiments pair as replicates 1 and 2 for the same time point experiments.

(4) Install JDK

(5) Install caper

## Step 3: Run the pipeline
The script submit.sh does the job. Output files are in /proj/yunligrp/users/tianyou/gRNA/NPC/source_file/ATAC-seq/output_slurm

## Step 4: summarize the results
Follow the instructions here to install croo: https://github.com/ENCODE-DCC/atac-seq-pipeline
Then run croo to get the reports.
Files generated are stored here: /proj/yunligrp/users/tianyou/gRNA/NPC/source_file/ATAC-seq/output_slurm/atac/03ff00ef-977b-4c1a-8e8f-407a6aa0fa56/align


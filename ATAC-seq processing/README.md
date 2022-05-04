# Processing raw ATAC-seq files for NPC
From the following paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6850896/

File location: 

NarrowPeak format: 	https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115046

SRA raw files: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA473806&f=assay_type_s%3An%3Aatac-seq%3Ac&o=acc_s%3Aa

Step 1: download the data

Select all and download the accession list. Then load the module sratoolkit on longleaf 
and follow the instructions at https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/ to configure SRA toolkit and download public data.
Following commands are used:

prefetch --option-file SRR_Acc_List.txt

for i in `seq 19 36`
do
echo SRR72304${i}
fasterq-dump --split-files SRR72304${i}.sra
done

Step 2: install the pipeline environment

Follow the instruction file to build the conda environment. The current anaconda version on longleaf could successfully build the environment.

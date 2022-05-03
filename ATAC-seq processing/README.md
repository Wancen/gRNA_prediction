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


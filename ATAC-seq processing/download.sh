#!/bin/bash

## prefetch data
prefetch --option-file SRR_Acc_List.txt

## convert to fastq files
for i in `seq 19 36`
do
echo SRR72304${i}
fasterq-dump --split-files SRR72304${i}.sra
done


#!/bin/bash
dir='/proj/yunligrp/users/tianyou/gRNA/codes/binary'

for fold in `seq 1 5`
do
for grp in pro enh
do
sbatch -N 1 -n 1 -p volta-gpu --mem=10g -t 03:00:00 --qos gpu_access --gres=gpu:1 --wrap="python -u ${dir}/CNN_binary_seq_topannot.py --fold ${fold} --grp ${grp}" -o ${grp}_seq_topannot_fold${fold}.txt
done
done

for fold in `seq 1 5`
do
for grp in pro enh
do
sbatch -N 1 -n 1 -p volta-gpu --mem=10g -t 03:00:00 --qos gpu_access --gres=gpu:1 --wrap="python -u ${dir}/CNN_binary_seq.py --fold ${fold} --grp ${grp}" -o ${grp}_seq_fold${fold}.txt
done
done

for fold in `seq 1 5`
do
for grp in pro enh
do
sbatch -N 1 -n 1 -p volta-gpu --mem=10g -t 03:00:00 --qos gpu_access --gres=gpu:1 --wrap="python -u ${dir}/CNN_binary_seqannot.py --fold ${fold} --grp ${grp}" -o ${grp}_seq_allannot_fold${fold}.txt
done
done

for fold in `seq 1 5`
do
for grp in pro enh
do
sbatch -N 1 -n 1 -p volta-gpu --mem=10g -t 03:00:00 --qos gpu_access --gres=gpu:1 --wrap="python -u ${dir}/CNN_binary_annot.py --fold ${fold} --grp ${grp}" -o ${grp}_annot_fold${fold}.txt
done
done

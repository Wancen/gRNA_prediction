#!/bin/bash
dir='/proj/yunligrp/users/tianyou/gRNA/codes/WT'

for fold in `seq 1 5`
do
sbatch -N 1 -n 1 -p volta-gpu --mem=10g -t 03:00:00 --qos gpu_access --gres=gpu:1 --wrap="python -u ${dir}/CNN_WT_continuous_seq.py --fold ${fold} --date Dec04" -o WT_seq_fold${fold}.txt
done

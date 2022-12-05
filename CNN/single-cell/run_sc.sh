#!/bin/bash

for i in `seq 1 5`
do
python -u CNN_sc_seqannot.py --fold $i --date Dec04 > sc_fold${i}.txt
done

for i in `seq 1 5`
do
python -u CNN_sc_seq.py --fold $i --date Dec04 > sc_seq_fold${i}.txt
done

for i in `seq 1 5`
do
python -u CNN_sc_annot.py --fold $i --date Dec04 > sc_annot_fold${i}.txt
done


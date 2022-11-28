#!/bin/bash

for i in `seq 1 5`
do
for cell in k562 ipsc npc
do
python -u mhc_binary_seqannot.py --fold $i --celltype ${cell} > ${cell}_fold${i}.txt
done
done

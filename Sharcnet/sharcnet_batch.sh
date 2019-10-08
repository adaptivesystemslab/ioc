#!/bin/bash

pwd

for g in `ls -d z_batch_iit*.sl`; do
  sqsubStr="sbatch $g"
  echo $sqsubStr
  $sqsubStr
  
  sleep 300s
done

#for f in `ls *.err`; do
#  rm $f
#done

#for f in `ls *.out`; do
#  rm $f
#done  

sleep 5s

squeue -u jf2lin
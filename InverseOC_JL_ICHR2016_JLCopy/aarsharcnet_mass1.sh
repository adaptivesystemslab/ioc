#!/bin/bash

dayCount=5d
mppCount=3G

cd script
pwd

for g in `ls -d sim_script3*`; do
  cd $g
  pwd
  
  for f in `ls sim*.m`; do
    sqsubStr="  sqsub -q serial -r $dayCount --mpp=$mppCount -o matlab_log_$f.log -i $f -j $f ../../../../uwmatlab"
    echo $sqsubStr
    $sqsubStr
  done
  
  for f in `ls *.m.log`; do
    rm $f
  done
  
  cd ..
done

sleep 1s
sqjobs -r -q
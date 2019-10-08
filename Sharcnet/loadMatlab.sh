#!/bin/bash
LD_PRELOAD=/usr/lib64/libstdc++.so.6
module load nixpkgs
module load matlab/2018b
chmod 775 ../../kalmanfilter/General_FKEKF/DynamicsModelMatlab/MatlabWrapper/DynamicsModelMatlab.mexa64
ldd -d ../../kalmanfilter/General_FKEKF/DynamicsModelMatlab/MatlabWrapper/DynamicsModelMatlab.mexa64

#matlab -nodisplay -nojvm

#!/bin/bash -l
#SBATCH --job-name=iit_test
#SBATCH --output=%x-%j.out
#SBATCH --error=%j.err
#SBATCH --account=def-dkulic   # adjust this to match the accounting group you are using to submit jobs
#SBATCH --time=1-00:00:00      # adjust this to match the walltime of your job
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # adjust this if you are using parallel commands
#SBATCH --mem=3000             # adjust this according to your the memory requirement per node you need
#SBATCH --mail-user=jf2lin@uwaterloo.ca # adjust this to match your email address
#SBATCH --mail-type=ALL

LD_PRELOAD=/usr/lib64/libstdc++.so.6
module load nixpkgs
module load matlab/2018b
chmod 775 /project/6001934/jf2lin/gitlab/kalmanfilter/General_FKEKF/DynamicsModelMatlab/MatlabWrapper/DynamicsModelMatlab.mexa64
ldd -d /project/6001934/jf2lin/gitlab/kalmanfilter/General_FKEKF/DynamicsModelMatlab/MatlabWrapper/DynamicsModelMatlab.mexa64

# without -sd, the starting dir is '/project/6001934/jf2lin'
srun matlab -nodisplay -nojvm -singleCompThread -r "IOCMain" -sd "/project/6001934/jf2lin/gitlab/expressive-ioc/InverseOC"

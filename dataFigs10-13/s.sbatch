#!/bin/bash
#
#SBATCH --job-name=agtest
#SBATCH --output=out.txt
#SBATCH --error=err.txt
#SBATCH --ntasks-per-node=28
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH -p gpu-long
#SBATCH --mail-type=ALL
#SBATCH --mail-user=saajid.chowdhury@stonybrook.edu

module load matlab/2023a

cd /gpfs/home/saachowdhury/agtest

matlab -batch "atomiongpu"

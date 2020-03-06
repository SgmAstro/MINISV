#!/bin/sh

#SBATCH --qos=regular
#SBATCH --time=1
#SBATCH --nodes=4
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell

source /project/projectdirs/desi/software/desi_environment.sh master

python run_rr_commands.py

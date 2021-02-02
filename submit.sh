#!/bin/bash
#SBATCH -J SM2  # Job name
#SBATCH -n 56  # Number of total cores
#SBATCH -N 1  # Number of nodes
#SBATCH --mem-per-cpu=2000# Memory pool for all cores in MB (see also --mem-per-cpu)
#SBATCH -o stdout # File to which STDOUT will be written %j is the job #
#SBATCH -p cpu
#SBATCH -A marom

module load cuda/9.0
module load intel/18.0.0.128
module load impi/2018
module load vasp/5.4.4_impi_cuda

python main.py

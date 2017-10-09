#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:00:30
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=3000
# Memory per node specification is in MB. It is optional.
# The default limit is 3000MB per core.
#SBATCH --job-name="snorreThreads"
#SBATCH --output=RPS.out
#SBATCH --mail-user=snorreob@ntnu.no
#SBATCH --mail-type=ALL

module load GCC
time ./serialRPS.out
time ./RPS_omp 8
time ./RPS_pthread 8

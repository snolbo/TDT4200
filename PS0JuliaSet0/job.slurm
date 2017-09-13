#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:00:30
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1200
# Memory per node specification is in MB. It is optional.
# The default limit is 3000MB per core.
#SBATCH --job-name="hello_test"
#SBATCH --output=test-srun.out
#SBATCH --mail-user=snorreob@ntnu.no
#SBATCH --mail-type=ALL


make
time ./julia_mpi 0.0 0.0

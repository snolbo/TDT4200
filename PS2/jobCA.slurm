#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:01:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --mem=3000
# Memory per node specification is in MB. It is optional.
# The default limit is 3000MB per core.
#SBATCH --job-name="JuliaSetMPISnolbo"
#SBATCH --output=RPS.out
#SBATCH --mail-user=snorreob@ntnu.no
#SBATCH --mail-type=ALL

module load GCC
module load OpenMPI
make
mpirun -np 4 ./RPSparallel

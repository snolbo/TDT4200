#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:00:30
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=500
# Memory per node specification is in MB. It is optional.
# The default limit is 3000MB per core.
#SBATCH --job-name="JuliaSetMPISnolbo"
#SBATCH --output=mpi_julia.out
#SBATCH --mail-user=snorreob@ntnu.no
#SBATCH --mail-type=ALL

module load GCC
module load OpenMPI
make
mpirun -np 4 ./julia_mpi -0.0 0.0

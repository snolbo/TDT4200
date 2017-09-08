#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:00:30
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=5
#SBATCH --mem=1200
# Memory per node specification is in MB. It is optional.
# The default limit is 3000MB per core.
#SBATCH --job-name="hello_test"
#SBATCH --output=test-srun.out
#SBATCH --mail-user=snorreob@ntnu.no
#SBATCH --mail-type=ALL


module load GCC
module load OpenMPI/2.0.2
make
time mpirun -n 10 ./collective

uname -a

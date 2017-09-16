#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=3000
# Memory per node specification is in MB. It is optional. 
# The default limit is 3000MB per core.
#SBATCH --job-name="hello_test"
#SBATCH --output=julia.out
#SBATCH --mail-user=snorreob@ntnu.no
#SBATCH --mail-type=ALL
 
 
make remake
./julia 0.258 0.01

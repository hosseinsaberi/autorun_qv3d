#!/bin/bash

# Standard output and error:
#SBATCH -o output-%j.out
#SBATCH -e std-%j.err

# Initial working directory:
#SBATCH -D ./

# Job Name:
#SBATCH -J qv3d

# Queue (Partition):
#SBATCH -p scarf
#SBATCH -n 320 

# Wall clock limit:
#SBATCH --time=01:00:00  # Example: 1 hour

# Run the program:
# ###############

module load FFTW.MPI/3.3.10-gompi-2023b
module load HDF5/1.14.3-gompi-2023b

srun -n  $SLURM_NTASKS ./qv3dMPIX.e v.ini

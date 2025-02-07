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
###SBATCH --time=01:00:00  # Example: 1 hour
#SBATCH --time=72:00:00  # Example: 3 days





# Load the modules:
# #################
module load FFTW.MPI/3.3.10-gompi-2023b
module load HDF5/1.14.3-gompi-2023b

# Record the start time
# #####################
echo "Job started at: $(date)" > job_times_${SLURM_JOB_ID}.txt

# Run QV3D job on clsuter
# #######################
srun -n  $SLURM_NTASKS ./qv3dMPIX.e v.ini

# Record the end time
# ###################
echo "Job ended at: $(date)" >> job_times-${SLURM_JOB_ID}.txt

# Calculate and record the run time
start_time=$(date -d "$(head -n 1 job_times_${SLURM_JOB_ID}.txt | cut -d ' ' -f 4-)" +%s)
end_time=$(date -d "$(tail -n 1 job_times_${SLURM_JOB_ID}.txt | cut -d ' ' -f 4-)" +%s)
run_time=$((end_time - start_time))
echo "Job run time: ${run_time} seconds" >> job_times-${SLURM_JOB_ID}.txt

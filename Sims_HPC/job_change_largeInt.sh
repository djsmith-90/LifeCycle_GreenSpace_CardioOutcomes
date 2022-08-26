#!/bin/bash    
                                                                
#SBATCH --job-name=change_largeInt
#SBATCH --partition=cpu
#SBATCH --nodes=1 # This is the number of nodes required for the job
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00 # (day-)hour:minute:second per job
#SBATCH --mem=5G
#SBATCH --array=1 # This is the number of datasets to run the script on
#SBATCH --error=/user/home/ds16565/LifeCycle/Error/Error_change_largeInt.txt # Error files (if any)
#SBATCH --output=/user/home/ds16565/LifeCycle/Output/Output_change_largeInt.txt # SLURM output file


# Change directory to save output in
cd "/user/home/ds16565/LifeCycle/Results/"

# Some other useful information
echo Running on hist "$(hostname)"
echo Time is "$(date)"
echo Directory is "$(pwd)"
echo Slurm job ID is "${SLURM_JOBID}"
echo This job runs on the following machines:
echo "${SLURM_JOB_NODELIST}"


# Load relevant software
module load languages/r/4.1.0

# Run the R script
srun Rscript /user/home/ds16565/LifeCycle/Scripts/change_largeInt.r


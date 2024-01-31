#!/bin/bash
#-------------------------------------------------------------------------------
#  BATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition hpc4   # for list of possibilities see http://docs.rnet.missouri.edu/Policy/partition
#SBATCH --nodes=1
#SBATCH --ntasks=1          # used for MPI codes, otherwise leave at '1'
#SBATCH --cpus-per-task=1  # cores per task
 
##SBATCH --exclusive        # using MPI with 90+% of the cores you should go exclusive
 
##SBATCH --mem=64G           # memory per node
#SBATCH --mem-per-cpu=1G   # memory per core (default is 1GB/core)
                            # mem and mem-per-cpu are mutually exclusive
 
#SBATCH --time 0-00:05:00     # days-hours:minutes
#SBATCH --qos=normal
#SBATCH --account=general   # investors will replace this with their account name
 
 
## labels and outputs
#SBATCH --job-name=compile_mrclean-greedy_job
#SBATCH --output=results-%j.out  # %j is the unique jobID
 
 
## notifications
#SBATCH --mail-user=kpsc59@mail.umsl.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL         # which type of notifications to send
#-------------------------------------------------------------------------------

# Load modules here:
module load boost

module list

# Compile program here
srun make

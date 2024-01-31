#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH --partition hpc4   # for list of possibilities see http://docs.rnet.missouri.edu/Policy/partition
#SBATCH --nodes=1
#SBATCH --ntasks=4          # used for MPI codes, otherwise leave at '1'
#SBATCH --cpus-per-task=1  # cores per task
 
##SBATCH --exclusive        # using MPI with 90+% of the cores you should go exclusive
 
##SBATCH --mem=64G           # memory per node
#SBATCH --mem-per-cpu=4G   # memory per core (default is 1GB/core)
                            # mem and mem-per-cpu are mutually exclusive
 
#SBATCH --time 0-05:00:00     # days-hours:minutes
#SBATCH --qos=normal
#SBATCH --account=climerlab   # investors will replace this with their account name
 
 
## labels and outputs
#SBATCH --job-name=run_sync_job
#SBATCH --output=results-%j.out  # %j is the unique jobID
 
 
## notifications
#SBATCH --mail-user=kpsc59@mail.umsl.edu  # email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL         # which type of notifications to send
#-------------------------------------------------------------------------------

# Load modules here:
export OMPI_MCA_btl_openib_if_include='qib0:1'
module load rss/rss-2020
module load openmpi/openmpi-3.1.3
module load boost

# Compile program here
srun /home/kpsc59/sync-greedy/sync sync.cfg 

#!/bin/bash
#SBATCH --job-name=ntm_battery
#SBATCH --nodes=1
#SBATCH --tasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --time=0-48:00:00
#SBATCH --array=0-6

module purge
module load intel/19.1.2
module load openmpi/intel/4.1.1

cd ../
cd $SLURM_ARRAY_TASK_ID

srun /scratch/mk8347/newlammps/src/lmp_mpi -in run_MD.lmp -screen screen.log

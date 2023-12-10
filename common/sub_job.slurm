#!/bin/bash

#SBATCH --job-name=ntm_battery
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10GB
#SBATCH --time=2-00:00:00

module purge
module load intel/19.1.2
module load openmpi/intel/4.1.1

srun lammps -in run_MD.lmp -screen screen.log

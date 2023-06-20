#!/bin/bash

#SBATCH --job-name=ntm_battery_2
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --time=0-06:00:00

module purge
module load intel/19.1.2
module load openmpi/intel/4.1.1

srun lammps -in run_MD.lmp -screen screen.log

### run your executable progra
######################

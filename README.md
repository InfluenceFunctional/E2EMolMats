
# LAMMPS_prep

A series of scripts for preparing and submitting LAMMPS simulations of various molecular solids.
Currently hardcoded for nicotinamide only, with options for vacancies, random reorientations and benzamide defects.


## 1. Installation
1. Requires following python packages: numpy, ase, scipy, ovito, pandas, mdanalysis. 
2. Requires also a local installation of moltemplate and LAMMPS compiled with the following functions:

|Name of Package|
|---------------|
|KSPACE         |
|MANYBODY       |
|MOLECULE       |
|OPT            |
|EXTRA-MOLECULE |
3. Edit paths in `generate_cluster_structures.py` and `main.py` to your specs.

## 2. Execution
1. Copy `configs/dev.py` to a new config.
2. Set desired simulation options. NOTE the script will loop over all possible combinations of list elements in the first half of the config, which can easily lead to a very large number of runs. See comments in `dev.py` for details. 
3. To run, import the desired config in `main.py` and run. NOTE: this will spawn new slurm jobs for each trajectory.  
4. Example script for slurm submission. NOTE: Module load calls may be different on your platform. 
```bash
#!/bin/bash

#SBATCH --job-name=lammps_md
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20GB
#SBATCH --time=0-04:00:00

module purge
module load intel/19.1.2
module load openmpi/intel/4.1.1

source activate your_python_environment

python main.py
```

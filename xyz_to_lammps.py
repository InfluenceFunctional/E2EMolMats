"""
generate an .xyz file and automatically prep LAMMPS inputs
"""

import warnings

warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')

import numpy as np
import os
from distutils.dir_util import copy_tree
from run_script import create_xyz_and_run_lammps

'''set head directory'''
head_dir = r'/home/mk8347/scratch/molecule_clusters/battery_11'
#head_dir = r'C:\Users\mikem\crystals\clusters\cluster_structures/battery_11'

crystals_path = r'/scratch/mk8347/molecule_clusters/CrystalStructures/'  #
#crystals_path = r'C:\Users\mikem\crystals\clusters\Leslie\CrystalStructures/'  #

if not os.path.exists(head_dir):
    os.mkdir(head_dir)
os.chdir(head_dir)
if not os.path.exists('common'):
    os.mkdir('common')
    copy_tree('../common', './common/')

cluster_sizes = [[4, 4, 4]]
temperatures = [200]
crystal_structures = ["NICOAM13"]
defect_rates = [0, 0.25, 0.5, 0.75, 1]

n_runs = len(cluster_sizes) * len(temperatures) * len(crystal_structures) * len(defect_rates)
run_nums = list(np.arange(1, n_runs + 1))

ind = 0
size_list = []
temp_list = []
crystal_list = []
defect_list = []
for i in range(len(cluster_sizes)):
    for j in range(len(temperatures)):
        for k in range(len(crystal_structures)):
            for l in range(len(defect_rates)):
                size_list.append(cluster_sizes[i])
                temp_list.append(temperatures[j])
                crystal_list.append(crystal_structures[k])
                defect_list.append(defect_rates[l])

for run_num, size, temp, crystal, defect in zip(run_nums, size_list, temp_list, crystal_list, defect_list):
    create_xyz_and_run_lammps(head_dir, run_num, crystals_path,
                              cluster_size=size,
                              print_steps=100,
                              run_time=int(1e5),
                              integrator='nosehoover',
                              box_type='p',
                              bulk_crystal=False,
                              min_inter_cluster_distance=500,
                              seed=1,
                              damping=str(100.0),
                              structure_identifier=crystal,
                              temperature=temp,
                              defect_rate=defect,
                              )

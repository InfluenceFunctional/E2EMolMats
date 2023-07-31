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
temperatures = [200, 300, 400]
crystal_structures = ["NICOAM13", "NICOAM17"]
gap_rates = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

n_runs = len(cluster_sizes) * len(temperatures) * len(crystal_structures) * len(gap_rates)
run_nums = list(np.arange(1, n_runs + 1))

ind = 0
size_list = []
temp_list = []
crystal_list = []
gap_list = []
for i in range(len(cluster_sizes)):
    for j in range(len(temperatures)):
        for k in range(len(crystal_structures)):
            for l in range(len(gap_rates)):
                size_list.append(cluster_sizes[i])
                temp_list.append(temperatures[j])
                crystal_list.append(crystal_structures[k])
                gap_list.append(gap_rates[l])

for run_num, size, temp, crystal, gap in zip(run_nums, size_list, temp_list, crystal_list, gap_list):
    create_xyz_and_run_lammps(head_dir, run_num, crystals_path,
                              cluster_size=size,
                              print_steps=1000,
                              run_time=int(1e7),
                              integrator='nosehoover',
                              box_type='s',
                              seed=1,
                              damping=str(100.0),
                              structure_identifier=crystal,
                              temperature=temp,
                              gap_rate=gap
                              )

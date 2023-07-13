"""
generate an .xyz file and automatically prep LAMMPS inputs
"""

import numpy as np
import os
from distutils.dir_util import copy_tree
from utils import create_xyz_and_run_lammps

'''set head directory'''
head_dir = r'/home/mk8347/scratch/molecule_clusters/bulk_reference'
#head_dir = r'C:\Users\mikem\crystals\clusters\cluster_structures/bulk_test'

crystals_path = r'/scratch/mk8347/molecule_clusters/CrystalStructures/'#
#crystals_path = r'C:\Users\mikem\crystals\clusters\Leslie\CrystalStructures/'  #

if not os.path.exists(head_dir):
    os.mkdir(head_dir)
os.chdir(head_dir)
if not os.path.exists('common'):
    os.mkdir('common')
    copy_tree('../common', './common/')

cluster_sizes = [[1, 1, 1]]
temperatures = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
crystal_structures = ["NICOAM13", "NICOAM17"]
n_runs = len(cluster_sizes) * len(temperatures) * len(crystal_structures)
run_nums = list(np.arange(1, n_runs + 1))

ind = 0
size_list = []
temp_list = []
crystal_list = []
for i in range(len(cluster_sizes)):
    for j in range(len(temperatures)):
        for k in range(len(crystal_structures)):
            size_list.append(cluster_sizes[i])
            temp_list.append(temperatures[j])
            crystal_list.append(crystal_structures[k])

for run_num, size, temp, crystal in zip(run_nums, size_list, temp_list, crystal_list):
    create_xyz_and_run_lammps(head_dir, run_num, crystals_path,
                              cluster_size=size,
                              print_steps=1000,
                              run_time=int(1e7),
                              integrator='nosehoover',
                              box_type='p',
                              seed=1,
                              damping=str(100.0),
                              structure_identifier=crystal,
                              temperature=temp,
                              )

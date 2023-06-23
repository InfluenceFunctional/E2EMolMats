"""
generate an .xyz file and automatically prep LAMMPS inputs
"""

import numpy as np
import os
from distutils.dir_util import copy_tree
from utils import create_xyz_and_run_lammps

'''set head directory'''
# head_dir = r'/home/mk8347/scratch/molecule_clusters/battery_7'
head_dir = r'C:\Users\mikem\crystals\clusters\cluster_structures/battery_7'

# crystals_path = r'/scratch/mk8347/molecule_clusters/CrystalStructures/'
crystals_path = r'C:\Users\mikem\crystals\clusters\Leslie\CrystalStructures/'  #

if not os.path.exists(head_dir):
    os.mkdir(head_dir)
os.chdir(head_dir)
if not os.path.exists('common'):
    os.mkdir('common')
    copy_tree('../common', './common/')

run_nums = list(np.arange(1, 19))
box_sizes = [50, 100, 150, 200, 400, 600] * 3
box_types = list(np.repeat(['p', 'f', 's'], 6))

for run_num, box_size, bound_type in zip(run_nums, box_sizes, box_types):
    create_xyz_and_run_lammps(head_dir, run_num, crystals_path,
                              min_inter_cluster_distance=box_size,
                              box_type=bound_type)

"""
generate an .xyz file and automatically prep LAMMPS inputs
"""

import numpy as np
import os
from distutils.dir_util import copy_tree
from utils import create_xyz_and_run_lammps

'''set head directory'''
head_dir = r'/home/mk8347/scratch/molecule_clusters/battery_7'
# head_dir = r'C:\Users\mikem\crystals\clusters\cluster_structures/battery_7'

crystals_path = r'/scratch/mk8347/molecule_clusters/CrystalStructures/'
# crystals_path = r'C:\Users\mikem\crystals\clusters\Leslie\CrystalStructures/'  #

if not os.path.exists(head_dir):
    os.mkdir(head_dir)
os.chdir(head_dir)
if not os.path.exists('common'):
    os.mkdir('common')
    copy_tree('../common', './common/')

run_nums = list(np.arange(1, 4))
integrator = 'langevin'
dampings = [str(1.0), str(10.0), str(100.0)]

for run_num, damp in zip(run_nums, dampings):
    create_xyz_and_run_lammps(head_dir, run_num, crystals_path,
                              integrator=integrator,
                              damping=damp)

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
head_dir = r'C:\Users\mikem\crystals\clusters\cluster_structures/benzamide_test2'

crystals_path = r'C:\Users\mikem\crystals\clusters\Leslie\CrystalStructures/'  #

if not os.path.exists(head_dir):
    os.mkdir(head_dir)
os.chdir(head_dir)
if not os.path.exists('common'):
    os.mkdir('common')
    copy_tree('../common', './common/')

create_xyz_and_run_lammps(head_dir, 1, crystals_path,
                          cluster_size=[4, 4, 4],
                          print_steps=1000,
                          run_time=int(1e7),
                          integrator='nosehoover',
                          box_type='p',
                          bulk_crystal=False,
                          min_inter_cluster_distance=500,
                          seed=1,
                          damping=str(100.0),
                          structure_identifier="NICOAM17",
                          temperature=200,
                          gap_rate=0,
                          defect_rate=0.5
                          )

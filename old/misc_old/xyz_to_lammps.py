"""
generate an .xyz file and automatically prep LAMMPS inputs
"""

import warnings

warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')

import numpy as np
import os
from distutils.dir_util import copy_tree
from e2emolmats.workflows.lammps_prepper import prep_lammps_inputs

'''set head directory'''
#head_dir = r'/home/mk8347/scratch/molecule_clusters/benzamide_test4'
head_dir = r'C:\Users\mikem\crystals\clusters\cluster_structures/misc_1'

#crystals_path = r'/scratch/mk8347/molecule_clusters/CrystalStructures/'  #
crystals_path = r'C:\Users\mikem\crystals\clusters\Leslie\CrystalStructures/'  #

if not os.path.exists(head_dir):
    os.mkdir(head_dir)
os.chdir(head_dir)
if not os.path.exists('md_data'):
    os.mkdir('md_data')
    copy_tree('../md_data', './md_data/')

cluster_sizes = [[3, 3, 3], [4, 4, 4], [5, 5, 5]]
temperatures = [300]
crystal_structures = ["NICOAM13", "NICOAM17"]
defect_rates = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5]

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
    prep_lammps_inputs(head_dir, run_num, crystals_path,
                       cluster_size=size,
                       print_steps=1000,
                       run_time=int(1e7),
                       integrator='nosehoover',
                       box_type='p',
                       bulk_crystal=False,
                       min_inter_cluster_distance=1000,
                       seed=1,
                       damping=str(100.0),
                       structure_identifier=crystal,
                       temperature=temp,
                       defect_rate=defect,
                       )

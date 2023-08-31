"""
generate an .xyz file and automatically prep LAMMPS inputs
"""

import itertools
import os
import sys
import warnings
from distutils.dir_util import copy_tree

warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')
from run_script import create_xyz_and_run_lammps

'''import run config'''
from configs.dev import batch_config

dynamic_configs = {key: value for key, value in batch_config.items() if isinstance(value, list)}
run_args = list(itertools.product(*list(dynamic_configs.values())))
dynamic_arg_keys = {key: i for i, key in enumerate(dynamic_configs.keys())}

n_runs = len(run_args)
print(f"Running {int(n_runs)} LAMMPS MD Jobs")

'''setup working directory'''
machine = batch_config['machine']  # 'cluster' or 'local'
if machine == 'local':
    #head_dir = r'/home/mkilgour/cluster_structures/' + batch_config['run_name']
    #crystals_path = r'/home/mkilgour/cluster_structures/CrystalStructures/'
    head_dir = r'C:\Users\mikem\crystals\clusters/cluster_structures/' + batch_config['run_name']
    crystals_path = r'C:\Users\mikem\crystals\clusters\Leslie/CrystalStructures/'
elif machine == 'cluster':
    head_dir = r'/home/mk8347/scratch/molecule_clusters/' + batch_config['run_name']
    crystals_path = r'/scratch/mk8347/molecule_clusters/CrystalStructures/'
else:
    print("Machine must be 'local' or 'cluster'")
    sys.exit()

if not os.path.exists(head_dir):
    os.mkdir(head_dir)

source_path = os.getcwd()
os.chdir(head_dir)

if not os.path.exists('common'):
    os.mkdir('common')
    copy_tree(source_path + '/common', './common/')  # copy from source

for run_num, run_config in enumerate(run_args):
    create_xyz_and_run_lammps(
        run_num=run_num,
        head_dir=head_dir,
        crystals_path=crystals_path,

        print_steps=batch_config['print_steps'],
        run_time=batch_config['run_time'],
        integrator=batch_config['integrator'],
        box_type=batch_config['box_type'],
        bulk_crystal=batch_config['bulk_crystal'],
        min_inter_cluster_distance=batch_config['min_inter_cluster_distance'],
        cluster_type=batch_config['cluster_type'],
        max_sphere_radius=batch_config['max_sphere_radius'],

        cluster_size=run_config[dynamic_arg_keys['cluster_size']],
        seed=run_config[dynamic_arg_keys['seed']],
        damping=run_config[dynamic_arg_keys['damping']],
        structure_identifier=run_config[dynamic_arg_keys['structure_identifier']],
        temperature=run_config[dynamic_arg_keys['temperature']],
        defect_rate=run_config[dynamic_arg_keys['defect_rate']],
        gap_rate=run_config[dynamic_arg_keys['gap_rate']],
        scramble_rate=run_config[dynamic_arg_keys['scramble_rate']],
    )

"""
generate an .xyz file and automatically prep LAMMPS inputs
"""

import warnings
import argparse
from utils import get_user_config, args2run_num, get_user_paths, generate_run_configs, \
    get_run_config, setup_workdir

warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')
from lammps_prepper import prep_LAMMPS_input

parser = argparse.ArgumentParser()
_, args = parser.parse_known_args()

'''Parse args, load config, get paths, setup workdir'''
run_num = args2run_num(args)
user_config = get_user_config(args)
batch_config = get_run_config(args)
run_args, dynamic_arg_keys = generate_run_configs(batch_config)
head_dir, crystals_path, ltemplify_path = get_user_paths(batch_config, user_config)
setup_workdir(head_dir)


def run_lammps_prepper(run_num, run_args, static_config, dynamic_arg_keys, ltemplify_path,head_dir, crystals_path):
    variable_config = run_args[run_num]
    run_config = static_config.copy()
    run_config.update({key: value for key, value in zip(dynamic_arg_keys.keys(), variable_config)})
    prep_LAMMPS_input(run_num, run_config,
                      ltemplify_path=ltemplify_path,
                      head_dir=head_dir,
                      crystals_path=crystals_path
                      )
'''generate structures and prep for MD'''
static_config = {key: value for key, value in batch_config.items() if not isinstance(value, list)}
if run_num is not None:  # for running in one at a time with run_num taken from args
    run_lammps_prepper(run_num, run_args, static_config, dynamic_arg_keys, ltemplify_path,head_dir, crystals_path)
else:  # for running sequentially in a loop
    for run_num, variable_config in enumerate(run_args):
        run_lammps_prepper(run_num, run_args, static_config, dynamic_arg_keys, ltemplify_path, crystals_path)

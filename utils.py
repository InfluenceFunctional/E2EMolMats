from argparse import Namespace
from generate_cluster_structures import generate_structure
from ovito.io import import_file, export_file
from template_scripts.initial_setup_for_ovito import initial_setup
from template_scripts.original_templify_to_runnable import templify_to_runnable
import ovito
from distutils.dir_util import copy_tree
import os

names_dict = {'1': 'H',  # rename for xyz export
              '8': 'H',
              '2': 'H',
              '6': 'N',
              '7': 'N',
              '4': 'C',
              '5': 'C',
              '3': 'O',
              }

ff_names_dict = {'1': 'ha',  # detailed atom types for analysis
                 '8': 'h4',
                 '2': 'hn',
                 '6': 'n',
                 '7': 'nb',
                 '4': 'c',
                 '5': 'ca',
                 '3': 'o',
                 }


def dict2namespace(data_dict: dict):
    """
    Recursively converts a dictionary and its internal dictionaries into an
    argparse.Namespace

    Parameters
    ----------
    data_dict : dict
        The input dictionary

    Return
    ------
    data_namespace : argparse.Namespace
        The output namespace
    """
    for k, v in data_dict.items():
        if isinstance(v, dict):
            data_dict[k] = dict2namespace(v)
        else:
            pass
    data_namespace = Namespace(**data_dict)

    return data_namespace


def create_xyz_and_run_lammps(head_dir, run_num, crystals_path,
                              cluster_type="supercell", structure_identifier="NICOAM13",
                              max_sphere_radius=20, cluster_size=None,
                              defect_rate=0, scramble_rate=0, gap_rate=0,
                              seed=1, min_inter_cluster_distance=200,
                              temperature=300, run_time=int(1e6),
                              print_steps=100, box_type='s',
                              integrator='langevin', damping: str = str(100.0)):
    # options
    if cluster_size is None:
        cluster_size = [3, 3, 3]

    workdir = head_dir + '/' + str(run_num)

    '''make new workdir'''
    if workdir is not None:
        if not os.path.exists(workdir):
            os.mkdir(workdir)
        os.chdir(workdir)
        '''copy in common elements'''
        copy_tree('../common', './')
    else:
        os.chdir(workdir)

    '''set temperature, run time, and print step in lmp file'''
    with open("run_MD.lmp") as f:
        newText = f.read().replace('_TEMP', str(temperature))
        newText = newText.replace('_RUNTIME', str(run_time))
        newText = newText.replace('_PRINTSTEPS', str(print_steps))
        newText = newText.replace('_SEED', str(seed))
        newText = newText.replace('_BOUND', str(box_type))
        newText = newText.replace('_DAMP', damping)
        if integrator == 'langevin':
            newText = newText.replace('#_LANGEVIN', '')
        elif integrator == 'nosehoover':
            newText = newText.replace('#_NOSE', '')

    with open("run_MD.lmp", "w") as f:
        f.write(newText)

    '''generate cluster structure'''
    xyz_filename = generate_structure(workdir, crystals_path, structure_identifier, cluster_type, max_sphere_radius,
                                      cluster_size, defect_rate, scramble_rate, gap_rate, seed, min_inter_cluster_distance)

    '''convert from .xyz to lammps datafile'''
    pipeline = import_file(xyz_filename)
    export_file(pipeline, '1.data', 'lammps/data', atom_style='full')

    '''prep for ovito bonds'''
    initial_setup(workdir, '1.data', '2.data')

    '''add bonds via ovito'''
    ovito.scene.load("nicotinamide_bond_session.ovito")
    pipeline = ovito.scene.pipelines[0]
    pipeline.source.load('2.data')
    export_file(pipeline, '3.data', 'lammps/data', atom_style='full')

    '''ltemplify'''
    os.system('ltemplify.py 3.data > 4.lt')  # .py on ltemplify required on cluster not windows

    '''make runnable'''
    templify_to_runnable(workdir, "4.lt", "3.data", "5.lt")

    '''run moltemplate and cleanup'''
    os.system("moltemplate.sh system.lt")
    os.system("cleanup_moltemplate.sh")

    '''optionally - directly run MD'''
    os.system("sbatch sub_job.slurm")

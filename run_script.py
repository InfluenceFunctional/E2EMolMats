import os
from distutils.dir_util import copy_tree
import numpy as np

import ovito
from ovito.io import import_file, export_file

from generate_cluster_structures import generate_structure
from template_scripts.initial_setup_for_ovito import initial_setup
from template_scripts.original_templify_to_runnable import templify_to_runnable
import subprocess


def create_xyz_and_run_lammps(head_dir, run_num, crystals_path, cluster_size,
                              cluster_type="supercell", structure_identifier="NICOAM13",
                              max_sphere_radius=None,
                              defect_rate=0, scramble_rate=0, gap_rate=0,
                              seed=1, min_inter_cluster_distance=500,
                              temperature=300, run_time=int(1e6),
                              print_steps=100, box_type='s', bulk_crystal=False,
                              integrator='langevin', damping: str = str(100.0)):
    """
    :param head_dir:
    :param run_num:
    :param crystals_path:
    :param cluster_size:
    :param cluster_type:
    :param structure_identifier:
    :param max_sphere_radius:
    :param defect_rate:
    :param scramble_rate:
    :param gap_rate:
    :param seed:
    :param min_inter_cluster_distance:
    :param temperature:
    :param run_time:
    :param print_steps:
    :param box_type:
    :param integrator:
    :param damping:
    :return:
    """

    '''make new workdir'''
    workdir = head_dir + '/' + str(run_num)
    if workdir is not None:
        if not os.path.exists(workdir):
            os.mkdir(workdir)
        os.chdir(workdir)
        '''copy in common elements'''
        copy_tree('../common', './')
    else:
        os.chdir(workdir)

    '''save run config'''
    run_config = {
        'head_dir': head_dir,
        'run_num': run_num,
        'crystals_path': crystals_path,
        'cluster_size': cluster_size,
        'cluster_type': cluster_type,
        'structure_identifier': structure_identifier,
        'max_sphere_radius': max_sphere_radius,
        'defect_rate': defect_rate,
        'scramble_rate': scramble_rate,
        'gap_rate': gap_rate,
        'seed': seed,
        'min_inter_cluster_distance': min_inter_cluster_distance,
        'temperature': temperature,
        'run_time': run_time,
        'print_steps': print_steps,
        'box_type': box_type,
        'bulk_crystal': bulk_crystal,
        'integrator': integrator,
        'damping': damping
    }
    np.save('run_config', run_config)

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
        elif integrator == 'npt':
            newText = newText.replace('#_NPT', '')

    with open("run_MD.lmp", "w") as f:
        f.write(newText)

    '''generate cluster structure'''
    xyz_filename = generate_structure(
        workdir, crystals_path, structure_identifier,
        cluster_type, max_sphere_radius,
        cluster_size, defect_rate, scramble_rate,
        gap_rate, seed, min_inter_cluster_distance,
        periodic_structure=bulk_crystal)

    '''convert from .xyz to lammps datafile'''
    pipeline = import_file(xyz_filename)
    export_file(pipeline, '1.data', 'lammps/data', atom_style='full')

    '''prep for ovito bonds'''
    initial_setup(workdir, '1.data', '2.data')

    '''add bonds via ovito'''
    ovito.scene.load("nicotinamide_bond_session_nico_ben_iso.ovito")
    pipeline = ovito.scene.pipelines[0]
    pipeline.source.load('2.data')
    export_file(pipeline, '3.data', 'lammps/data', atom_style='full')

    '''ltemplify'''
    if r'Users\mikem' in workdir:  # if we are Mike's local windows machine (which
        ltemplify_path = r"C:\Users\mikem\miniconda3\envs\LAMMPS_runs\lib\site-packages\moltemplate\ltemplify.py"
    else:  # works on linux
        ltemplify_path = subprocess.getoutput("which ltemplify.py")

    os.system(f'python {ltemplify_path} 3.data > 4.lt')  # .py on ltemplify required on cluster not windows

    '''make runnable'''
    templify_to_runnable(workdir, "4.lt", "3.data", "5.lt")

    '''run moltemplate and cleanup'''
    os.system("moltemplate.sh system.lt")
    os.system("cleanup_moltemplate.sh")

    # '''optionally - directly run MD'''
    #os.system("sbatch sub_job.slurm")

"""
generate an .xyz file and automatically prep LAMMPS inputs
"""

from generate_cluster_structures import generate_structure
from ovito.io import import_file, export_file
from template_scripts.initial_setup_for_ovito import initial_setup
from template_scripts.original_templify_to_runnable import templify_to_runnable
import ovito
import os
from distutils.dir_util import copy_tree

'''set head directory'''
head_dir = r'/home/mk8347/scratch/molecule_clusters/battery_6'
# head_dir = r'C:\Users\mikem\crystals\clusters\cluster_structures/battery_5'

crystals_path = r'/scratch/mk8347/molecule_clusters/CrystalStructures/'
# crystals_path = r'C:\Users\mikem\crystals\clusters\Leslie\CrystalStructures/'  #

if not os.path.exists(head_dir):
    os.mkdir(head_dir)
os.chdir(head_dir)
if not os.path.exists('common'):
    os.mkdir('common')
    copy_tree('../common', './common/')

run_nums = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
gaps = [0, 0.05, 0.1, 0.15, 0.2, 0, 0.05, 0.1, 0.15, 0.2]
structure_ids = ["NICOAM13", "NICOAM13", "NICOAM13", "NICOAM13", "NICOAM13",
                 "NICOAM17", "NICOAM17", "NICOAM17", "NICOAM17", "NICOAM17"]

for run_num, gap, structure_id in zip(run_nums, gaps, structure_ids):
    # options
    workdir = head_dir + '/' + str(run_num)
    structure_identifier = structure_id
    cluster_type = "supercell"  # "supercell" or "spherical"
    max_sphere_radius = 20
    cluster_size = [3, 3, 3]  # size of supercell to be generated [a,b,c]
    defect_rate = 0  # fraction of molecules which will be switched from nicotinamide to benzamide
    scramble_rate = 0  # fraction of molecules with scrambled orientations
    gap_rate = gap  # fraction of molecules deleted
    seed = 1  # must be positive integer
    min_inter_cluster_distance = 200  # angstroms
    temperature = 300
    run_time = 10000000
    print_steps = 1000

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

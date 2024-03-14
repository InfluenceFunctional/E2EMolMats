import os
from distutils.dir_util import copy_tree
import numpy as np

import ovito
from ovito.io import import_file, export_file
from ovito.data import *
from ovito.pipeline import *
from ovito.modifiers import *
from argparse import Namespace

from generate_cluster_structures import generate_structure
from template_scripts.initial_setup_for_ovito import initial_setup as initial_setup_nicotinamide
from template_scripts.acridine_initial_for_ovito import initial_setup as initial_setup_acridine
from template_scripts.original_templify_to_runnable import templify_to_runnable as templify_to_runnable_nicotinamide
from template_scripts.acridine_original_templify_to_runnable import \
    templify_to_runnable as templify_to_runnable_acridine
from template_scripts.moltemp_final import moltemp_final as moltemp_final_nicotinamide
from template_scripts.acridine_moltemp_final_ver2 import moltemp_final as moltemp_final_acridine
import subprocess


def settings_final(change_atom_style=True):
    file = 'system.in.init'
    data = open(file, 'r')
    New_data = open(f'new_{file}', 'w')

    lines = data.readlines()
    if change_atom_style:
        lines[0] = "atom_style full2\n"
        #    for i in range(8):
        #        line = lines[i]
        #        #split_line = line.split('/')
        #        #split_line[3] = split_line[3].replace('charmm', 'long')
        #        #new_line = '/'.join(split_line)
        #        New_data.write(new_line)
        #
        #    New_data.write('pair_coeff 9 9 lj/charmm/coul/long 0.0860 3.39966950842\npair_coeff 10 10 lj/charmm/coul/long 0.0860 3.39966950842\n')

    for i in range(0, len(lines)):
        New_data.write(lines[i])

    New_data.close()


def settings_final_nico(change_atom_style=True):
    file = 'system.in.init'
    data = open(file, 'r')
    New_data = open(f'new_{file}', 'w')

    lines = data.readlines()
    if change_atom_style:
        lines[0] = "atom_style full\n"
        #    for i in range(8):
        #        line = lines[i]
        #        #split_line = line.split('/')
        #        #split_line[3] = split_line[3].replace('charmm', 'long')
        #        #new_line = '/'.join(split_line)
        #        New_data.write(new_line)
        #
        #    New_data.write('pair_coeff 9 9 lj/charmm/coul/long 0.0860 3.39966950842\npair_coeff 10 10 lj/charmm/coul/long 0.0860 3.39966950842\n')

    for i in range(0, len(lines)):
        New_data.write(lines[i])

    New_data.close()


def create_xyz_and_run_lammps(run_config):
    """
    main working script
    """
    config = Namespace(**run_config)

    '''make new workdir'''
    workdir = config.head_dir + '/' + str(config.run_num)
    if os.path.exists(workdir):
        pass
    else:
        if not os.path.exists(workdir):
            os.mkdir(workdir)
        os.chdir(workdir)
        '''copy in common elements'''
        copy_tree('../common', './')

        '''save run config'''

        np.save('run_config', run_config)

        print("============================")
        print("Generating Structure")
        print("============================")
        '''generate cluster structure'''
        xyz_filename, melt_inds = generate_structure(
            workdir,
            config.crystals_path,
            config.structure_identifier,
            config.cluster_type,
            config.max_sphere_radius,
            config.cluster_size,
            config.defect_rate,
            config.defect_type,
            config.scramble_rate,
            config.gap_rate,
            config.seed,
            config.min_inter_cluster_distance,
            config.min_lattice_length,
            periodic_structure=config.bulk_crystal,
            prep_crystal_in_melt=config.prep_crystal_in_melt)

        '''set temperature, run time, and print step in lmp file'''
        with (open("run_MD.lmp") as f):
            newText = f.read()

            if config.integrator.lower() == 'langevin':
                newText = newText.replace('#_LANGEVIN', '')
            elif config.integrator.lower() == 'nosehoover':
                newText = newText.replace('#_NOSE', '')
            elif config.integrator.lower() == 'npt':
                newText = newText.replace('#_NPT', '')
            if config.bulk_crystal:
                newText = newText.replace('#_KSPACE', '')
            if config.prep_crystal_in_melt:
                newText = newText.replace('#_MELT_PREP', '')
                newText = newText.replace('_EQUIL_TIME', str(config.equil_time))
                newText = newText.replace('_MELT_TEMP', str(config.melt_temperature))
                newText = newText.replace('_MELT_START_IND', str(melt_inds.melt_start_ind))
                newText = newText.replace('_MELT_END_IND', str(melt_inds.melt_end_ind))
                newText = newText.replace('_CRYSTAL_START_IND', str(melt_inds.crystal_start_ind))
                newText = newText.replace('_CRYSTAL_END_IND', str(melt_inds.crystal_end_ind))

            newText = newText.replace('_TEMP_SAMPLE', str(config.temperature))
            newText = newText.replace('_RUNTIME', str(config.run_time))
            newText = newText.replace('_PRINTSTEPS', str(config.print_steps))
            newText = newText.replace('_SEED', str(config.seed))
            newText = newText.replace('_BOUND', str(config.box_type))
            newText = newText.replace('_DAMP', config.damping)

            if 'nicotinamide' in config.structure_identifier:
                newText = newText.replace('#_NICOTINAMIDE', '')

            elif 'acridine' in config.structure_identifier:
                newText = newText.replace('#_ACRIDINE', '')

            if config.ramp_temperature:
                newText = newText.replace('_INIT_TEMP', str(int(1)))
            else:
                newText = newText.replace('_INIT_TEMP', str(config.temperature))

        with open("run_MD.lmp", "w") as f:
            f.write(newText)

        print("============================")
        print("Converting to lammps and ovito analysis")
        print("============================")

        '''convert from .xyz to lammps datafile'''
        pipeline = import_file(xyz_filename)
        export_file(pipeline, '1.data', 'lammps/data', atom_style='full')

        '''prep for ovito bonds'''
        if 'nicotinamide' in config.structure_identifier:
            pipeline.source.load('1.data')
            create_bonds_modifier = CreateBondsModifier(mode=CreateBondsModifier.Mode.VdWRadius)
            pipeline.modifiers.append(create_bonds_modifier)

            # CreateBondsModifier.Mode.VdWRadius
            # pipeline.modifiers.append(CreateBondsModifier.Mode.VdWRadius())
            # pipeline.modifiers.append(CreateBondsModifier())
            export_file(pipeline, '2.data', 'lammps/data', atom_style='full')
            initial_setup_nicotinamide(workdir, '2.data', '3.data')
            # ovito.scene.load("nicotinamide_bond_session_nico_ben_iso.ovito")
            # create_bonds_modifier = CreateBondsModifier(mode=CreateBondsModifier.Mode.VdWRadius)
            # pipeline.modifiers.append(create_bonds_modifier)
            # export_file(pipeline, '3.data', 'lammps/data', atom_style='full')

        elif 'acridine' in config.structure_identifier:
            initial_setup_acridine(workdir, '1.data', '2.data')
            pipeline.source.load('2.data')
            pipeline.modifiers.append(CreateBondsModifier(cutoff=1.7))
            export_file(pipeline, '3.data', 'lammps/data', atom_style='full')
            # ovito.scene.load("acridine_ovito.ovito")
            # pipeline = import_file("2.data")
            # pipeline.modifiers.append(CreateBondsModifier(cutoff = 1.7))
            # pipeline.create_bond(a, b, type=None, pbcvec=None)

        '''add bonds via ovito'''
        #        pipeline = ovito.scene.pipelines[0]
        # pipeline.source.load('2.data')
        #        pipeline.modifiers.append(CreateBondsModifier(cutoff = 1.7))
        # modbonds = CreateBondsModifier(mode = CreateBondsModifier.Mode.VdWRadius)
        # pipeline.modifiers.append(modbonds)
        #        pipeline.modifiers.append(CreateBondsModifier.Mode.VdWRadius)
        # pipeline.modifiers.append(CreateBondsModifier.Mode.VdWRadius)
        # pipeline.modifiers.append(CreateBondsModifier.Mode.VdWRadius)
        #        pipeline.modifiers.append(CreateBondsModifier.Mode.VdWRadius)

        print("============================")
        print("Ltemplifying")
        print("============================")

        '''ltemplify'''
        if r'Users\mikem' in workdir:  # if we are Mike's local windows machine (which
            ltemplify_path = r"C:\Users\mikem\miniconda3\envs\LAMMPS_runs\lib\site-packages\moltemplate\ltemplify.py"
        else:  # works on linux
            ltemplify_path = subprocess.getoutput("unset -f which; which ltemplify.py")
        #            print(ltemplify_path)
        #        print("###################################################################################")
        #        print(ltemplify_path)
        # os.system(f'{ltemplify_path} 3.data > 4.lt')  # .py on ltemplify required on cluster not windows
        os.system("~/.local/bin/ltemplify.py 3.data > 4.lt")  # .py on ltemplify required on cluster not windows

        print("============================")
        print("Templify to runnable")
        print("============================")

        '''make runnable'''
        if 'nicotinamide' in config.structure_identifier:
            templify_to_runnable_nicotinamide(workdir, "4.lt", "3.data", "5.lt")
        elif 'acridine' in config.structure_identifier:
            templify_to_runnable_acridine(workdir, "4.lt", "3.data", "5.lt")

        print("============================")
        print("Running Moltemplate")
        print("============================")

        '''run moltemplate and cleanup'''
        os.system("~/.local/bin/moltemplate.sh system.lt")

        print("============================")
        print("Moltemplate cleanup")
        print("============================")

        os.system("~/.local/bin/cleanup_moltemplate.sh")

        print("============================")
        print("Indexing cleanup")
        print("============================")

        if 'nicotinamide' in config.structure_identifier:  # these are in fact un-used in nicotinamide runs
            moltemp_final_nicotinamide(workdir)  # Daisuke final indexing cleanup
            settings_final_nico(True)  # adjust pairs to be Daisuke-friendly

        elif 'acridine' in config.structure_identifier:
            moltemp_final_acridine(workdir)  # Daisuke final indexing cleanup
            settings_final(True)  # adjust pairs to be Daisuke-friendly

        print("============================")
        print("Submitting LAMMPS run")
        print("============================\n")

        # '''optionally - directly run MD''' # NOTE this no longer works on cluster due to Singularity issue - must use the batch_sub_lmp.sh in common to submit after all templates are built
        # os.system("/opt/slurm/bin/sbatch sub_job.slurm")

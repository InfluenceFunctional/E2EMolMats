import os
from distutils.dir_util import copy_tree
import numpy as np

from ovito.io import import_file, export_file
from ovito.modifiers import *
from argparse import Namespace

from generate_cluster_structures import generate_structure
from template_scripts.initial_setup import initial_setup
from template_scripts.initial_setup_for_ovito import initial_setup as initial_setup_nicotinamide
from template_scripts.acridine_initial_for_ovito import initial_setup as initial_setup_acridine
from template_scripts.moltemp_final import moltemp_final
from template_scripts.original_templify_to_runnable import templify_to_runnable as templify_to_runnable_nicotinamide
from template_scripts.acridine_original_templify_to_runnable import \
    templify_to_runnable as templify_to_runnable_acridine

from template_scripts.utils import update_atom_style_in_settings, generate_MD_script


def prep_LAMMPS_input(run_num, config_i, ltemplify_path, head_dir, crystals_path):
    """
    main working script
    """
    config = Namespace(**config_i)
    molecule_name = config.structure_identifier.split('/')[0]  # molecule name is 1st half of structure identifier with format "molecule/form"

    '''make new workdir'''
    workdir = head_dir + '/' + str(run_num)
    if os.path.exists(workdir):
        print(f"{workdir} already exists - skipping run")
        pass
    else:
        os.mkdir(workdir)
        os.chdir(workdir)
        '''copy in common elements'''
        copy_tree('../common', './')
        np.save('run_config', config)

        print("============================")
        print("Generating Structure")
        print("============================")
        '''generate cluster structure'''
        xyz_filename, melt_inds, molind2name = generate_structure(
            crystals_path,
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
        generate_MD_script(config, melt_inds)

        print("============================")
        print("Converting to lammps and ovito analysis")
        print("============================")

        '''convert from .xyz to lammps datafile'''
        pipeline = import_file(xyz_filename)
        export_file(pipeline, '1.data', 'lammps/data', atom_style='full')

        '''prep atom types and make bonds'''
        initial_setup('1.data', '2.data', molecule_name, molind2name)

        pipeline.source.load('2.data')
        create_bonds_modifier = CreateBondsModifier(mode=CreateBondsModifier.Mode.VdWRadius)
        pipeline.modifiers.append(create_bonds_modifier)
        export_file(pipeline, '2.data', 'lammps/data', atom_style='full')
        initial_setup_nicotinamide(workdir, '2.data', '3.data')
        #
        # initial_setup_acridine(workdir, '1.data', '2.data')
        # pipeline.source.load('2.data')
        # pipeline.modifiers.append(CreateBondsModifier(cutoff=1.7))
        # export_file(pipeline, '3.data', 'lammps/data', atom_style='full')

        print("============================")
        print("Ltemplifying")
        print("============================")

        '''ltemplify'''
        #ltemplify_path = subprocess.getoutput("unset -f which; which ltemplify.py") # alternate method
        os.system(f"{ltemplify_path} 3.data > 4.lt")  # .py on ltemplify required on cluster not windows

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

        moltemp_final(workdir, config.atom_style, molind2name)  # final indexing cleanup
        update_atom_style_in_settings(atom_style=config.atom_style)

        if config.submit_lammps_slurm:
            print("============================")
            print("Submitting LAMMPS run")
            print("============================\n")

            # '''optionally - directly run MD''' # will not work from within a Singularity instance
            # use instead batch_sub_lmp.sh in /common to submit after all templates are built
            os.system("/opt/slurm/bin/sbatch sub_job.slurm")



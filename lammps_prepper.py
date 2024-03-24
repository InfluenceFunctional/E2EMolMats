import os
from distutils.dir_util import copy_tree
import numpy as np

from ovito.io import import_file, export_file
from ovito.modifiers import *
from argparse import Namespace

from generate_cluster_structures import generate_structure
from template_scripts.initial_setup import initial_setup
from template_scripts.moltemp_final import moltemp_final
from template_scripts.templify_to_runnable import templify_to_runnable

from template_scripts.utils import update_atom_style_in_settings, generate_MD_script


def prep_lammps_inputs(run_num, config_i, ltemplify_path, head_dir, crystals_path):
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
        if head_dir.split('/')[-1] == 'dev':  # index the next highest integer dev directory
            ind = 0
            found_new_workdir = False
            while not found_new_workdir:
                if os.path.exists(workdir + "_" + str(ind)):
                    ind += 1
                else:
                    found_new_workdir = True

        os.mkdir(workdir)
        os.chdir(workdir)
        '''copy in common elements'''
        copy_tree('../common', './')

        print("============================")
        print("Generating Structure")
        print("============================")
        '''generate cluster structure'''
        xyz_filename, melt_inds, molind2name = generate_structure(  # TODO save molind2name in run_config
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

        config_dict = config.__dict__
        config_dict.update({'molind2name_dict': molind2name})
        np.save('run_config', config_dict)

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
        create_bonds_modifier = CreateBondsModifier(cutoff=1.7)
        #create_bonds_modifier = CreateBondsModifier(mode=CreateBondsModifier.Mode.VdWRadius)
        pipeline.modifiers.append(create_bonds_modifier)
        export_file(pipeline, '3.data', 'lammps/data', atom_style='full')

        print("============================")
        print("Ltemplifying")
        print("============================")

        '''ltemplify'''
        #ltemplify_path = subprocess.getoutput("unset -f which; which ltemplify.py") # alternate method
        os.system(f"{ltemplify_path} 3.data > 4.lt")

        print("============================")
        print("Templify to runnable")
        print("============================")

        '''make runnable'''
        templify_to_runnable('4.lt', '3.data', '5.lt',
                             molecule_name)

        print("============================")
        print("Running Moltemplate")
        print("============================")

        '''run moltemplate and cleanup'''  # todo change to user-config path
        os.system("~/.local/bin/moltemplate.sh system.lt -nocheck")  # nocheck means it will skip over missing @bond type issues

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



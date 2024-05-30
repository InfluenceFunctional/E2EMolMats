import os
from distutils.dir_util import copy_tree
import numpy as np

from ovito.io import import_file, export_file
from ovito.modifiers import *
from argparse import Namespace

from e2emolmats.common.generate_cluster_structures import generate_structure
from e2emolmats.datafile_processing.initial_setup import atom_type_renumbering
from e2emolmats.datafile_processing.moltemp_final import moltemp_final
from e2emolmats.datafile_processing.templify_to_runnable import templify_to_runnable
from e2emolmats.datafile_processing.utils import update_atom_style_in_settings, generate_MD_script
from e2emolmats.md_data.constants import MOLECULE_NUM_ATOMS, ATOM_TYPES, MOLECULE_SYM_INDICES


def prep_lammps_inputs(run_num, config_i, ltemplify_path, head_dir, crystals_path):
    """
    main working script
    """
    config = Namespace(**config_i)
    # molecule name is 1st half of structure identifier with format "molecule/form"
    molecule_name = config.structure_identifier.split('/')[0]

    '''make new workdir'''
    if head_dir.split('/')[-1] != 'dev':  # if not dev
        workdir = head_dir + '/' + str(run_num)
        if os.path.exists(workdir):  # if it's been done, skip run
            print(f"{workdir} already exists - skipping run")
            return False
    else:  # if dev, find next empty slot
        workdir = head_dir + '/'
        ind = 0
        found_new_workdir = False
        while not found_new_workdir:
            if os.path.exists(workdir + str(ind)):
                ind += 1
            else:
                found_new_workdir = True
                workdir += str(ind)
    if not os.path.exists(workdir):
        os.mkdir(workdir)
    os.chdir(workdir)
    '''copy in md_data elements'''
    copy_tree('../md_data', './')

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
        prep_crystal_in_melt=config.prep_crystal_in_melt,
        prep_melt_interface=config.prep_melt_interface,
        prep_bulk_melt=config.prep_bulk_melt,
        melt_interface_direction=config.pressure_direction)

    config_dict = config.__dict__
    config_dict.update({'molind2name_dict': molind2name})
    config_dict.update({'molecule_num_atoms_dict': MOLECULE_NUM_ATOMS})
    config_dict.update({'molecule_atom_types': ATOM_TYPES})
    config_dict.update({'molecule_sym_indices': MOLECULE_SYM_INDICES})
    config_dict.update({'melt_indices': melt_inds})
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
    atom_type_renumbering('1.data', '2.data',
                          molecule_name, molind2name, config.defect_rate, config.defect_type)

    pipeline.source.load('2.data')
    create_bonds_modifier = CreateBondsModifier(cutoff=1.7, intra_molecule_only=True, prevent_hh_bonds=True)

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
    os.system(
        "~/.local/bin/moltemplate.sh system.lt -nocheck")  # nocheck means it will skip over missing @bond type issues

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
        # use instead batch_sub_lmp.sh in /md_data to submit after all templates are built
        os.system("/opt/slurm/bin/sbatch sub_job.slurm")

import numpy as np

from e2emolmats.md_data.constants import MOLECULE_NUM_ATOM_TYPES


def update_atom_style_in_settings(atom_style='full'):
    file = 'system.in.init'
    data = open(file, 'r')
    New_data = open(f'new_{file}', 'w')

    lines = data.readlines()
    lines[-1] = f"atom_style {atom_style}\n"

    for i in range(0, len(lines)):
        New_data.write(lines[i])

    New_data.close()


def generate_MD_script(config, melt_inds):
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
        if config.prep_crystal_in_melt or config.prep_melt_interface:
            newText = newText.replace('#_MELT_PREP', '')
            newText = newText.replace('#_CRYSTAL_PREP', '')
            newText = newText.replace('_EQUIL_TIME', str(config.equil_time))
            newText = newText.replace('_MELT_TEMP', str(config.melt_temperature))
            newText = newText.replace('_MELT_START_IND', str(melt_inds.melt_start_ind))
            newText = newText.replace('_MELT_END_IND', str(melt_inds.melt_end_ind))
            newText = newText.replace('_CRYSTAL_START_IND', str(melt_inds.crystal_start_ind))
            newText = newText.replace('_CRYSTAL_END_IND', str(melt_inds.crystal_end_ind))
        elif config.prep_bulk_melt:  # exclusive with above
            newText = newText.replace('#_MELT_PREP', '')
            newText = newText.replace('_EQUIL_TIME', str(config.equil_time))
            newText = newText.replace('_MELT_TEMP', str(config.melt_temperature))
            newText = newText.replace('_MELT_START_IND', str(melt_inds.melt_start_ind))  # melt inds are simply the full cluster
            newText = newText.replace('_MELT_END_IND', str(melt_inds.melt_end_ind))

        newText = newText.replace('_TEMP_SAMPLE', str(config.temperature))
        newText = newText.replace('_RUNTIME', str(config.run_time))
        newText = newText.replace('_PRINTSTEPS', str(config.print_steps))
        newText = newText.replace('_SEED', str(config.seed))
        newText = newText.replace('_BOUND', str(config.box_type))
        newText = newText.replace('_DAMP', config.damping)

        # if 'nicotinamide' in config.structure_identifier:
        #     if config.atom_style == 'full2':
        #         newText = newText.replace('#_NICOTINAMIDE_sym', '')
        #     else:
        #         newText = newText.replace('#_NICOTINAMIDE_no_sym', '')
        #
        # elif 'acridine' in config.structure_identifier:
        #     if config.atom_style == 'full2':
        #         newText = newText.replace('#_ACRIDINE_sym', '')
        #     else:
        #         newText = newText.replace('#_ACRIDINE_no_sym', '')

        if config.ramp_temperature:
            newText = newText.replace('_INIT_TEMP', str(config.init_temperature))
            newText = newText.replace('#_EQUIL_BEFORE_RAMP', '')
            newText = newText.replace('_EQUIL_TIME', str(config.equil_time))
        else:
            newText = newText.replace('_INIT_TEMP', str(config.temperature))

        # dump modify
        molecule = config.structure_identifier.split('/')[0]
        num_atom_types = MOLECULE_NUM_ATOM_TYPES[molecule]
        type_string = ''
        for ind in range(num_atom_types):
            type_string = type_string + f' {ind + 1}'
        newText = newText.replace('#_DUMP_MODIFY', '  dump_modify       d2 element' + type_string)

    with open("run_MD.lmp", "w") as f:
        f.write(newText)

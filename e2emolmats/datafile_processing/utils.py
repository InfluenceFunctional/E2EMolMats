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
        new_text = f.read()

        if config.integrator.lower() == 'langevin':
            new_text = new_text.replace('#_LANGEVIN', '')
        elif config.integrator.lower() == 'nosehoover':
            new_text = new_text.replace('#_NOSE', '')
        elif config.integrator.lower() == 'npt':
            new_text = new_text.replace('#_NPT', '')
            new_text = new_text.replace('#_p_dim', str(config.pressure_direction))
        if config.bulk_crystal:
            new_text = new_text.replace('#_KSPACE', '')
        if config.prep_crystal_in_melt or config.prep_melt_interface:
            new_text = new_text.replace('#_MELT_PREP', '')
            new_text = new_text.replace('#_CRYSTAL_IN_MELT_PREP', '')
            new_text = new_text.replace('_EQUIL_TIME', str(config.equil_time))
            new_text = new_text.replace('_MELT_TEMP', str(config.melt_temperature))
            new_text = new_text.replace('_MELT_START_IND', str(melt_inds.melt_start_ind))
            new_text = new_text.replace('_MELT_END_IND', str(melt_inds.melt_end_ind))
            new_text = new_text.replace('_CRYSTAL_START_IND', str(melt_inds.crystal_start_ind))
            new_text = new_text.replace('_CRYSTAL_END_IND', str(melt_inds.crystal_end_ind))
        elif config.prep_bulk_melt:  # exclusive with above
            new_text = new_text.replace('#_MELT_PREP', '')
            new_text = new_text.replace('_EQUIL_TIME', str(config.equil_time))
            new_text = new_text.replace('_MELT_TEMP', str(config.melt_temperature))
            new_text = new_text.replace('_MELT_START_IND', str(melt_inds.melt_start_ind))  # melt inds are simply the full cluster
            new_text = new_text.replace('_MELT_END_IND', str(melt_inds.melt_end_ind))

        new_text = new_text.replace('_TEMP_SAMPLE', str(config.temperature))
        new_text = new_text.replace('_RUNTIME', str(config.run_time))
        new_text = new_text.replace('_PRINTSTEPS', str(config.print_steps))
        new_text = new_text.replace('_SEED', str(config.seed))
        new_text = new_text.replace('_BOUND', str(config.box_type))
        new_text = new_text.replace('_DAMP', config.damping)

        # if 'nicotinamide' in config.structure_identifier:
        #     if config.atom_style == 'full2':
        #         new_text = new_text.replace('#_NICOTINAMIDE_sym', '')
        #     else:
        #         new_text = new_text.replace('#_NICOTINAMIDE_no_sym', '')
        #
        # elif 'acridine' in config.structure_identifier:
        #     if config.atom_style == 'full2':
        #         new_text = new_text.replace('#_ACRIDINE_sym', '')
        #     else:
        #         new_text = new_text.replace('#_ACRIDINE_no_sym', '')

        if config.ramp_temperature:
            new_text = new_text.replace('_INIT_TEMP', str(config.init_temperature))
            new_text = new_text.replace('#_EQUIL_BEFORE_RAMP', '')
            new_text = new_text.replace('_EQUIL_TIME', str(config.equil_time))
        else:
            new_text = new_text.replace('_INIT_TEMP', str(config.temperature))

        # dump modify
        molecule = config.structure_identifier.split('/')[0]
        if config.defect_rate > 0:
            num_atom_types = MOLECULE_NUM_ATOM_TYPES[molecule + '+' + config.defect_type]
        else:
            num_atom_types = MOLECULE_NUM_ATOM_TYPES[molecule]
        if molecule == 'acridine':
            num_atom_types = 4  # todo this is hardcoded for new GRACE FF

        type_string = ''
        for ind in range(num_atom_types):
            type_string = type_string + f' {ind + 1}'
        new_text = new_text.replace('#_DUMP_MODIFY', '  dump_modify       d2 element' + type_string)

    with open("run_MD.lmp", "w") as f:
        f.write(new_text)

import os
from e2emolmats.md_data.constants import MOLECULE_SYM_INDICES


def moltemp_final(workdir: str, atom_style: str, molind2name: dict) -> None:
    os.chdir(workdir)
    original_filename = "system.data"
    new_filename = "new_" + original_filename
    original_data = open(original_filename, 'r')
    new_data = open(new_filename, 'w')

    # first adjust atom_tyle
    line = original_data.readline()
    while line != "Atoms  # full\n":
        new_data.write(line)
        line = original_data.readline()
    new_data.write(f"Atoms  # {atom_style}\n")

    # adjust atom type legend  # todo check with Daisuke if we need to do this
    line = original_data.readline()
    atom_counter = 1

    if atom_style == 'full2':
        # full2 includes symmetry info
        while line != "\n":
            line = line.split(' ')

            # get type of current molecule
            mol_ind = line[1]
            sym_index_dict = MOLECULE_SYM_INDICES[molind2name[mol_ind]]

            # add symmetry tag
            line[3] = ' ' + str(sym_index_dict[atom_counter])
            newline = ' '.join(line) + '\n'
            new_data.write(newline)

            atom_counter += 1
            if atom_counter == len(sym_index_dict):  # this only works if molecules are processed in blocks
                atom_counter = 1
                mol_ind += 1

    # write the rest of the file
    new_data.write(line)
    while line:
        line = original_data.readline()
        new_data.write(line)

    # finish up
    new_data.close()
    original_data.close()

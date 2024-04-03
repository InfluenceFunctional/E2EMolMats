from e2emolmats.md_data.constants import (MOLECULE_ATOM_TYPES_MASSES,
                                          MOLECULE_NUM_ATOM_TYPES, MOLECULE_NUM_ATOMS,
                                          ATOM_TYPES, ATOM_CHARGES)


def initial_setup(original_filename, new_filename, molecule_name, molind2name):
    datafile = open(original_filename, 'r')
    new_datafile = open(new_filename, 'w')
    lines = datafile.readlines()

    num_atom_types = MOLECULE_NUM_ATOM_TYPES[molecule_name]

    mol_counter = 0
    atom_counter = 1
    hit_masses = False
    hit_atoms = False
    printed_masses = False
    printed_atoms = False
    new_mol = True
    inside_atoms_block = False
    end_of_atom = False
    element_index_dict = {}
    for counter, line in enumerate(lines):
        # get element indices from ovito export
        if hit_masses:
            if '# ' in line:
                splitline = line.split(' ')
                element_symbol = splitline[-1]
                if len(element_symbol) <= 2:
                    element_index_dict[element_symbol] = splitline[0]  # index assigned to this element

        if 'Masses' in line:
            hit_masses = True
        if 'Atoms' in line:
            hit_atoms = True

        if counter == 3:  # line 4 in the file, counter indexes from 0
            new_datafile.write(f"{num_atom_types} atom types\n")

        elif hit_masses and not hit_atoms and not printed_masses:
            new_datafile.write(MOLECULE_ATOM_TYPES_MASSES[molecule_name])
            new_datafile.write("\n")
            printed_masses = True

        elif hit_atoms and not printed_atoms:
            new_datafile.write(line)  # prints "Atoms  # full"
            new_datafile.write("\n")
            printed_atoms = True

        elif printed_atoms and line != '\n' and not end_of_atom:
            inside_atoms_block = True
            line = line.split()
            if new_mol:  # get information to print the given molecule
                mol_counter += 1
                molecule_type = molind2name[mol_counter]
                mol_num_atoms = MOLECULE_NUM_ATOMS[molecule_type]
                atom_counter = 1
                new_mol = False
                
            line[1] = str(mol_counter)
            if atom_counter == mol_num_atoms:  # reset and check for new molecule type on the next line
                new_mol = True
            
            write_molecule(line, atom_counter, new_datafile, molecule_type)
            atom_counter += 1

        elif inside_atoms_block and line == '\n':
            inside_atoms_block = False
            new_datafile.write(line)
            end_of_atom = True

        elif not hit_masses:  # do nothing
            new_datafile.write(line)

    new_datafile.close()
    datafile.close()


def write_molecule(line, atom_counter, new_datafile, molecule_type):
    line[2] = str(ATOM_TYPES[molecule_type][atom_counter])
    line[3] = str(ATOM_CHARGES[molecule_type][atom_counter])
    new_datafile.write(' '.join(line) + "\n")

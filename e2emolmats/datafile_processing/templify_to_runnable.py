from e2emolmats.md_data.constants import MOLECULE_SHORTHAND, FF_PATH_DICT


def templify_to_runnable(original_lt_path, original_data_path, new_lt_path, molecule_name):
    '''
    ff_paths gaff2_nicotinamid_long.lt, gaff2_acridine.lt
    '''
    molecule_shorthand = MOLECULE_SHORTHAND[molecule_name]
    ff_name = FF_PATH_DICT[molecule_name]

    data = open(original_lt_path, 'r')
    new_lt = open(new_lt_path, 'w')
    lt_lines = data.readlines()

    data_file = open(original_data_path, 'r')
    lt_file = open("system.lt", 'w')
    lt_file.write("import \"" + new_lt_path + "\"  # <- defines the molecule type.\n\n\n")
    lt_file.write("# Periodic boundary conditions:\nwrite_once(\"Data Boundary\") {\n")
    [data_file.readline() for _ in range(6)]

    data_file_line = data_file.readline()
    lt_file.write("   " + data_file_line)
    data_file_line = data_file.readline()
    lt_file.write("   " + data_file_line)
    '''xyz box dimensions'''
    data_file_line = data_file.readline()
    lt_file.write("   " + data_file_line)
    data_file_line = data_file.readline()
    lt_file.write("   " + data_file_line)
    data_file_line = data_file.readline()
    lt_file.write("   " + data_file_line)

    data_file_line = data_file.readline()
    if 'xy xz yz' in data_file_line:  # if there is a non-orthogonal component, write it
        lt_file.write("   " + data_file_line)

    lt_file.write("}\n\n")
    lt_file.write(f"# Create \"{molecule_name.capitalize()}\" molecules\n# rotated and moved to give polymorph {molecule_shorthand} = new {molecule_name.capitalize()}\n") # TODO do we need this?
    lt_file.write(f"{molecule_shorthand} = new {molecule_name.capitalize()}")

    counter = 0
    end_data_atom = False
    data_bond_checker = 0
    end_data_bond = 0
    hit_write_atoms_block = False
    for ind, line in enumerate(lt_lines):
        counter = counter + 1
        if counter < 5:
            new_lt.write(line)
        elif counter == 5:
            new_lt.write(line)
            new_lt.write(f"import \"{ff_name}\"\n{molecule_name.capitalize()}" + " inherits GAFF2 {\n")
        elif counter < 10:
            new_lt.write(line)
        elif counter >= 10 and not hit_write_atoms_block:
            if 'write("Data Atoms") {' in line:
                hit_write_atoms_block = True
                new_lt.write(line)
            else:
                new_lt.write(line)
        elif not end_data_atom:
            line = line.split()
            if line[0] == "}":
                new_lt.write("}")
                end_data_atom = True
            else:
                new_lt.write("  " + line[0] + " " + line[1] + " " + line[2] + " " + line[3] + " " + line[4] + " " + line[5] + " " + line[6] + "\n")
        elif data_bond_checker == 0:
            if line == "write(\"Data Bonds\") {\n":
                data_bond_checker = 1
                new_lt.write("write(\"Data Bond List\") {\n")
            else:
                new_lt.write(line)
        elif end_data_bond == 0:
            line = line.split()
            if line[0] == "}":
                new_lt.write("}\n")
                end_data_bond = 1
            else:
                new_lt.write("  " + line[0] + " " + line[2] + " " + line[3] + "\n")
    new_lt.write("}")

    new_lt.close()
    data.close()
    lt_file.close()


#templify_to_runnable(original_lt_path, original_data_path, new_lt_path, molecule_name)

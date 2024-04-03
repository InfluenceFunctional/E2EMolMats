# import os
#
# os.chdir(r'C:\\Users\\mikem\\crystals\\clusters\\cluster_structures\\benzamide_test2\\3')


def residue_cleanup():
    f = open('system.data', "r")
    text = f.read()
    lines = text.split('\n')
    lines_rstrip = [line.rstrip("\n") for line in lines]
    f.close()

    new_file = open('fixed_system.data', 'w')

    mol_type_ind = []
    atom_counter = 0
    mol_counter = 0
    new_mol = True
    counter = 0
    skip = True
    last_line_skip = False

    for ind, i in enumerate(lines):
        counter = counter + 1
        if skip:
            new_file.write(i)
            new_file.write("\n")

        if not skip:
            line = i.split()
            if line == []:
                skip = True
                last_line_skip = False
                new_file.write("\n")

        if not skip:
            line = i.split()
            if new_mol:  # If I am reading a new molecule, I will go here
                mol_lines_3rd = lines_rstrip[counter + 1].split()  # checking the 3rd atom
                mol_lines_4th = lines_rstrip[counter + 2].split()  # checking the 4th atom
                mol_counter += 1
                if mol_lines_3rd[2] == '7':  # ring nitrogen in Benzamide
                    mol_flag = 0
                    mol_atom_num = 15
                elif mol_lines_3rd[2] == '5' and mol_lines_4th[2] == '5':  # ring carbons
                    mol_flag = 1
                    mol_atom_num = 16
                elif mol_lines_3rd[2] == '5' and mol_lines_4th[2] == '7':  # shifted ring nitrogen
                    mol_flag = 2
                    mol_atom_num = 15
                # else:  # regularly triggers erroneously
                #     print(" mol_lines_3rd[2] = " + mol_lines_3rd[2])
                #     print("Something in the atom info is wrong. It is not nicotinamide, benzamide, or isonicotinamide")
                #     aa = 1
                #     exit()

                atm_counter = 1
                new_mol = False
            if atm_counter == mol_atom_num:  # reset and check for new molecule type
                new_mol = True

            line[1] = str(mol_counter)
            new_line = ' '.join(line)
            new_file.write(new_line + "\n")

            atm_counter += 1

        if last_line_skip:
            skip = False

        if skip and i == "Atoms  # full":
            last_line_skip = True

    new_file.close()

    aa = 0

# In this file, you convert data file with coordinate information
import os


def initial_setup(workdir, original, New):
    os.chdir(workdir)
    data = open(original, 'r')
    New_data = open(New, 'w')
    lines = data.readlines()

    lines_rstrip = [line.rstrip("\n") for line in lines]

    counter = 0
    mol_counter = 0
    atm_counter = 1
    hit_masses = False
    hit_atoms = False
    printed_masses = False
    printed_atoms = False
    mol_atom_num = 23  # for acridine
    new_mol = True
    for ind, i in enumerate(lines):
        counter = counter + 1
        if 'Masses' in i:
            hit_masses = True
        if 'Atoms' in i:
            hit_atoms = True

        if counter == 4:
            New_data.write("3 atom types\n")
        elif not hit_masses:
            New_data.write(i)
        elif hit_masses and not hit_atoms and not printed_masses:
            # elif counter == 11:
            # New_data.write("Masses\n\n1 1.008  # ha\n2 1.008  # h4\n3 1.008  # hn\n4 14.01  # n\n5 14.01  # nb\n6 12.01  # c\n7 12.01  # ca\n8 12.01  # ca1\n9 12.01  # ca2\n10 16.0  # o\n")
            # New_data.write("Masses\n\n1 1.008  # ha\n2 1.008  # h4\n3 1.008  # hn\n4 14.01  # n\n5 14.01  # nb\n6 12.01  # c\n7 12.01  # ca\n8 16.00  # o\n")
            New_data.write("Masses\n\n1 14.01  # nb\n2 12.01  # ca\n3 1.008  # ha\n")
            New_data.write("\n")
            # New_data.write(i)
            printed_masses = True
        elif hit_atoms and not printed_atoms:
            # elif counter == 17:
            New_data.write(i)
            New_data.write("\n")
            printed_atoms = True
        elif printed_atoms and i != '\n':
            # elif counter >= 19:
            line = i.split()
            if new_mol:  # mol_checker != line[1]:  # If I am reading a new molecule, I will go here
                mol_lines_3rd = lines_rstrip[counter + 1].split()  # checking the 3rd atom
                mol_lines_4th = lines_rstrip[counter + 2].split()  # checking the 4th atom
                mol_checker = mol_lines_3rd[1]
                mol_counter += 1
                atm_counter = 1
                new_mol = False
            line[1] = mol_counter
            if atm_counter == mol_atom_num:  # reset and check for new molecule type
                new_mol = True

            write_acridine(line, atm_counter, New_data)

            atm_counter += 1

    New_data.close()
    data.close()


def write_acridine(line, atm_counter, New_data):
    if atm_counter == 1:
        line[2] = 1
        line[3] = -0.693636
    elif atm_counter == 2:
        line[2] = 2
        line[3] = 0.534823
    elif atm_counter == 3:
        line[2] = 2
        line[3] = -0.002672
    elif atm_counter == 4:
        line[2] = 2
        line[3] = -0.20686
    elif atm_counter == 5:
        line[2] = 2
        line[3] = -0.002672
    elif (atm_counter == 6):
        line[2] = 2
        line[3] = 0.534823
    elif (atm_counter == 7):
        line[2] = 2
        line[3] = -0.274751
    elif (atm_counter == 8):
        line[2] = 2
        line[3] = -0.100851
    elif (atm_counter == 9):
        line[2] = 2
        line[3] = -0.140261
    elif (atm_counter == 10):
        line[2] = 2
        line[3] = -0.187088
    elif (atm_counter == 11):
        line[2] = 2
        line[3] = -0.187088
    elif (atm_counter == 12):
        line[2] = 2
        line[3] = -0.140261
    elif (atm_counter == 13):
        line[2] = 2
        line[3] = -0.100851
    elif (atm_counter == 14):
        line[2] = 2
        line[3] = -0.274751
    elif (atm_counter == 15):
        line[2] = 3
        line[3] = 0.151995
    elif (atm_counter == 16):
        line[2] = 3
        line[3] = 0.152509
    elif (atm_counter == 17):
        line[2] = 3
        line[3] = 0.128637
    elif (atm_counter == 18):
        line[2] = 3
        line[3] = 0.130168
    elif (atm_counter == 19):
        line[2] = 3
        line[3] = 0.133735
    elif (atm_counter == 20):
        line[2] = 3
        line[3] = 0.133735
    elif (atm_counter == 21):
        line[2] = 3
        line[3] = 0.130168
    elif (atm_counter == 22):
        line[2] = 3
        line[3] = 0.128637
    elif (atm_counter == 23):
        line[2] = 3
        line[3] = 0.152509

    for j in range(0, 7):
        New_data.write(str(line[j]) + ' ')
    New_data.write("\n")
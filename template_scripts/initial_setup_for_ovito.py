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
    mol_flag = -1  # 0 = nicotinamide, 1 = benzamide, 2 = isonicotinamide.
    hit_masses = False
    hit_atoms = False
    printed_masses = False
    printed_atoms = False
    n_benzamide_atoms = 16
    n_nicotinamide_atoms = 15
    new_mol = True
    oxy = -1;
    nit = -1;
    car = -1;
    for ind, i in enumerate(lines):
        counter = counter + 1
        if "# C" in i:
            car = i[0];
        if "# O" in i:
            oxy = i[0];
        if "# N" in i:
            nit = i[0];
        if 'Masses' in i:
            hit_masses = True
        if 'Atoms' in i:
            hit_atoms = True

        if counter == 4:
            New_data.write("8 atom types\n")
        elif not hit_masses:
            New_data.write(i)
        elif hit_masses and not hit_atoms and not printed_masses:
            # elif counter == 11:
            #New_data.write("Masses\n\n1 1.008  # ha\n2 1.008  # h4\n3 1.008  # hn\n4 14.01  # n\n5 14.01  # nb\n6 12.01  # c\n7 12.01  # ca\n8 12.01  # ca1\n9 12.01  # ca2\n10 16.0  # o\n")
            New_data.write("Masses\n\n1 1.008  # ha\n2 1.008  # h4\n3 1.008  # hn\n4 14.01  # n\n5 14.01  # nb\n6 12.01  # c\n7 12.01  # ca\n8 16.00  # o\n")
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
            if new_mol: #mol_checker != line[1]:  # If I am reading a new molecule, I will go here
                mol_lines_3rd = lines_rstrip[counter + 1].split()  # checking the 3rd atom
                mol_lines_4th = lines_rstrip[counter + 2].split()  # checking the 4th atom
                mol_checker = mol_lines_3rd[1]
                mol_counter += 1
                if mol_lines_3rd[2] == nit:  # ring nitrogen in Benzamide
                    mol_flag = 0
                    mol_atom_num = n_nicotinamide_atoms
                elif mol_lines_3rd[2] == car and mol_lines_4th[2] == car:  # ring carbons
                    mol_flag = 1
                    mol_atom_num = n_benzamide_atoms
                elif mol_lines_3rd[2] == car and mol_lines_4th[2] == nit:  # shifted ring nitrogen
                    mol_flag = 2
                    mol_atom_num = n_nicotinamide_atoms
                else:
                    print(" mol_lines_3rd[2] = " + mol_lines_3rd[2])
                    print("Something in the atom info is wrong. It is not nicotinamide, benzamide, or isonicotinamide")
                    aa = 1
                    exit()

                atm_counter = 1
                new_mol = False
            line[1] = mol_counter;
            if atm_counter == mol_atom_num:  # reset and check for new molecule type
                new_mol = True

            if mol_flag == 0:  # case of nicotinamide
                write_nicotinamide(line, atm_counter, New_data)

            if mol_flag == 1:  # case of benzamide
                write_benzamide(line, atm_counter, New_data)

            if mol_flag == 2:  # case of isonicotinamide
                write_isonicotinamide(line, atm_counter, New_data)

            atm_counter += 1

    New_data.close()
    data.close()


def write_nicotinamide(line, atm_counter, New_data):
    if (atm_counter == 1):
        line[2] = 7
        line[3] = -0.326615
        
    elif (atm_counter == 2):
        line[2] = 7
        line[3] = 0.380568
        
    elif (atm_counter == 3):
        line[2] = 5
        line[3] = -0.585364
        
    elif (atm_counter == 4):
        line[2] = 7
        line[3] = 0.384166
        
    elif (atm_counter == 5):
        line[2] = 7
        line[3] = -0.38538
        
    elif (atm_counter == 6):
        line[2] = 7
        line[3] = 0.173788
        
    elif (atm_counter == 7):
        line[2] = 6
        line[3] = 0.6054
        
    elif (atm_counter == 8):
        line[2] = 8
        line[3] = -0.479634
        
    elif (atm_counter == 9):
        line[2] = 4
        line[3] = -0.779885
        
    elif (atm_counter == 10):
        line[2] = 3
        line[3] = 0.357505
        
    elif (atm_counter == 11):
        line[2] = 3
        line[3] = 0.357505
        
    elif (atm_counter == 12):
        line[2] = 2
        line[3] = 0.027993
        
    elif (atm_counter == 13):
        line[2] = 2
        line[3] = 0.034858
        
    elif (atm_counter == 14):
        line[2] = 1
        line[3] = 0.157212
        
    elif (atm_counter == 15):
        line[2] = 1
        line[3] = 0.077882
        
    for j in range(0, 7):
        New_data.write(str(line[j]) + ' ')
    New_data.write("\n")


def write_benzamide(line, atm_counter, New_data):
    if (atm_counter == 1):
        line[2] = 7
        line[3] = -0.073516
        
    elif (atm_counter == 2):
        line[2] = 7
        line[3] = -0.081906
        
    elif (atm_counter == 3):
        line[2] = 7
        line[3] = -0.123099
        
    elif (atm_counter == 4):
        line[2] = 7  # Change the name of ca2 to ca2/n
        line[3] = -0.104813
        
    elif (atm_counter == 5):
        line[2] = 7
        line[3] = -0.123099
        
    elif (atm_counter == 6):
        line[2] = 7
        line[3] = -0.081906
        
    elif (atm_counter == 7):
        line[2] = 6
        line[3] = 0.58899
        
    elif (atm_counter == 8):
        line[2] = 8
        line[3] = -0.47625
        
    elif (atm_counter == 9):
        line[2] = 4
        line[3] = -0.787316
        
    elif (atm_counter == 10):
        line[2] = 3
        line[3] = 0.350274
        
    elif (atm_counter == 11):
        line[2] = 3
        line[3] = 0.350274
        
    elif (atm_counter == 12):
        line[2] = 2
        line[3] = 0.102542
        
    elif (atm_counter == 13):
        line[2] = 2
        line[3] = 0.119669
        
    elif (atm_counter == 14):
        line[2] = 1
        line[3] = 0.117946
        
    elif (atm_counter == 15):
        line[2] = 1
        line[3] = 0.119669
        
    elif (atm_counter == 16):
        line[2] = 1
        line[3] = 0.102542
        
    for j in range(0, 7):
        New_data.write(str(line[j]) + ' ')
    New_data.write("\n")


def write_isonicotinamide(line, atm_counter, New_data):
    if (atm_counter == 1):
        line[2] = 7
        line[3] = 0.228047
        
    elif (atm_counter == 2):
        line[2] = 7
        line[3] = -0.429772
        
    elif (atm_counter == 3):
        line[2] = 7
        line[3] = 0.429792
        
    elif (atm_counter == 4):
        line[2] = 7
        line[3] = -0.608913
        
    elif (atm_counter == 5):
        line[2] = 7
        line[3] = 0.429792
        
    elif (atm_counter == 6):
        line[2] = 7
        line[3] = -0.429772
        
    elif (atm_counter == 7):
        line[2] = 6
        line[3] = 0.560881
        
    elif (atm_counter == 8):
        line[2] = 8
        line[3] = -0.462231
        
    elif (atm_counter == 9):
        line[2] = 4
        line[3] = -0.809865
        
    elif (atm_counter == 10):
        line[2] = 3
        line[3] = 0.366217
        
    elif (atm_counter == 11):
        line[2] = 3
        line[3] = 0.366217
        
    elif atm_counter == 12:
        line[2] = 2
        line[3] = 0.153507
        
    elif atm_counter == 13:
        line[2] = 2
        line[3] = 0.026297
        
    elif (atm_counter == 14):
        line[2] = 1
        line[3] = 0.026297
        
    elif (atm_counter == 15):
        line[2] = 1
        line[3] = 0.153507
        
    for j in range(0, 7):
        New_data.write(str(line[j]) + ' ')
    New_data.write("\n")

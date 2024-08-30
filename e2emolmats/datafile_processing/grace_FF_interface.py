import os
from pathlib import Path
from shutil import copyfile


def grace_FF_preprocessing(input_filename, output_filename):
    grace_dir = Path('acridine_GRACE')
    ###########################
    # Hard code for acridine
    natm = 23  # number of atoms in a molecule
    atom_type = 4
    nbond = 25
    nangle = 40
    ndi = 138
    nimp = 0
    grace_atom_type = [1, 1, 1, 2, 3, 2, 1, 1, 1, 1, 2, 1, 2, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4]

    ###########################

    # initialization
    bond_type = 0
    angle_type = 0
    di_type = 0
    improper_type = 0
    nmol = 0
    lmp_ind = []
    gra_ind = []
    bond_lst = []
    angle_lst = []
    di_lst = []
    imp_lst = []

    # Get the charge info (w/ GRACE index)
    charge_f = open(os.path.join(grace_dir, "charge.txt"), "r")
    line = charge_f.readline()
    GRACEcharge = line.split()
    charge_f.close()

    #Get Acridine index used in prep and grace
    with open(os.path.join(grace_dir, "Acr_lam_to_grace_index.txt"), "r") as f:
        for line in f:
            txt = line.split()
            lmp_ind.append(int(txt[0]))
            gra_ind.append(int(txt[1]))
    #with open("bond_test.txt","r") as f:
    with open(os.path.join(grace_dir, "bond.txt"), "r") as f:
        for line in f:  # num, bond type, atom1, atom2
            bond_lst.append(line)
    with open(os.path.join(grace_dir, "angle.txt"), "r") as f:
        for line in f:
            angle_lst.append(line)
    with open(os.path.join(grace_dir, "di.txt"), "r") as f:
        for line in f:
            di_lst.append(line)
    with open(os.path.join(grace_dir, "imp.txt"), "r") as f:
        for line in f:
            imp_lst.append(line)

    # Read settings file. Need these for data file. Also trying to make it less hard coded. You're welcome future me!
    with open(os.path.join(grace_dir, "system.in.settings"), "r") as f:
        for line in f:
            txt = line.split()
            if txt[0] == "bond_coeff":
                bond_type += 1
            elif txt[0] == "angle_coeff":
                angle_type += 1
            elif txt[0] == "dihedral_coeff":
                di_type += 1
            elif txt[0] == "improper_coeff":
                improper_type += 1

    #old_data = open("new_acridine_single.data","r");
    ##################################################Please review file location of this part
    old_data = open(input_filename, "r")
    new_data = open(output_filename, "w")
    #Slip junk
    old_data.readline()
    old_data.readline()
    line = old_data.readline()
    txt = line.split()
    natom = int(txt[0])
    nmol = int(natom / natm)
    #skip old atom_type
    old_data.readline()

    new_data.write("# LAMMPS data file written by MK and DK\n\n")
    new_data.write(str(natom) + " atoms\n")
    new_data.write(str(int(nbond * nmol)) + " bonds\n")
    new_data.write(str(int(nangle * nmol)) + " angles\n")
    new_data.write(str(int(ndi * nmol)) + " dihedrals\n")
    new_data.write(str(int(nimp * nmol)) + " impropers\n\n")
    new_data.write(str(int(atom_type)) + " atom types\n")
    new_data.write(str(int(bond_type)) + " bond types\n")
    new_data.write(str(int(angle_type)) + " angle types\n")
    new_data.write(str(int(di_type)) + " dihedral types\n")
    new_data.write(str(int(improper_type)) + " improper types\n\n")

    while not (len(txt) == 4 and txt[-1] == "xhi"):
        line = old_data.readline()
        txt = line.split()
    new_data.write(line)

    while not (len(txt) == 1 and txt[0] == "Masses"):
        line = old_data.readline()
        txt = line.split()
        new_data.write(line)
    new_data.write("\n1 12.0107\n2 12.0107\n3 14.0067\n4 1.00794\n\n")
    while not (len(txt) >= 1 and txt[0] == "Atoms"):
        line = old_data.readline()
        txt = line.split()
    #Write Atom
    new_data.write(line)
    line = old_data.readline()
    #Write space
    new_data.write(line)
    #First line of atom geo
    line = old_data.readline()
    txt = line.split()
    #Read Atom part. Need to change charge info
    while (len(txt) != 0):
        this_lmp_index = int(txt[0]) % 23 - 1
        this_grace_index = gra_ind[this_lmp_index]
        #    print(this_lmp_index)
        #    print(this_grace_index)
        #    print(gra_ind)
        #    exit()
        txt[3] = GRACEcharge[this_grace_index]
        txt[2] = str(grace_atom_type[this_grace_index])
        new_line = ' '.join(txt)
        new_data.write(new_line)
        new_data.write("\n")
        line = old_data.readline()
        txt = line.split()

    #Make Bonds parts
    if (nbond > 0): new_data.write("\nBonds\n\n")
    bond_count = 1
    for i in range(0, nbond):
        #for i in range(0,1):
        txt = bond_lst[i].split()
        new_tixt = txt[:]
        #if you have more than one molecule, use this loop to create other bonds
        for j in range(0, nmol):
            new_txt = txt[:]
            #put appropriate i and j atom for each bond
            grace_index_i = int(txt[2]) % (natm + 1) - 1
            grace_index_j = int(txt[3]) % (natm + 1) - 1
            lmp_index_i = gra_ind.index(grace_index_i)
            lmp_index_j = gra_ind.index(grace_index_j)
            new_txt[0] = str(bond_count)
            bond_count += 1
            new_txt[2] = str((lmp_index_i + 1) + natm * j)
            new_txt[3] = str((lmp_index_j + 1) + natm * j)
            new_data.write(" ".join(new_txt))
            new_data.write("\n")

    if (nangle > 0): new_data.write("\nAngles\n\n")
    angle_count = 1
    for i in range(0, nangle):
        txt = angle_lst[i].split()
        new_tixt = txt[:]
        #if you have more than one molecule, use this loop to create other angles
        for j in range(0, nmol):
            new_txt = txt[:]
            #put appropriate i and j atom for each angle
            grace_index_i = int(txt[2]) % (natm + 1) - 1
            grace_index_j = int(txt[3]) % (natm + 1) - 1
            grace_index_k = int(txt[4]) % (natm + 1) - 1
            lmp_index_i = gra_ind.index(grace_index_i)
            lmp_index_j = gra_ind.index(grace_index_j)
            lmp_index_k = gra_ind.index(grace_index_k)
            new_txt[0] = str(angle_count)
            angle_count += 1
            new_txt[2] = str((lmp_index_i + 1) + natm * j)
            new_txt[3] = str((lmp_index_j + 1) + natm * j)
            new_txt[4] = str((lmp_index_k + 1) + natm * j)
            new_data.write(" ".join(new_txt))
            new_data.write("\n")

    if (ndi > 0): new_data.write("\nDihedrals\n\n")
    di_count = 1
    for i in range(0, ndi):
        txt = di_lst[i].split()
        new_tixt = txt[:]
        # if you have more than one molecule, use this loop to create other dihedrals
        for j in range(0, nmol):
            new_txt = txt[:]
            # put appropriate i and j atom for each dihedrals
            grace_index_i = int(txt[2]) % (natm + 1) - 1
            grace_index_j = int(txt[3]) % (natm + 1) - 1
            grace_index_k = int(txt[4]) % (natm + 1) - 1
            grace_index_l = int(txt[5]) % (natm + 1) - 1
            lmp_index_i = gra_ind.index(grace_index_i)
            lmp_index_j = gra_ind.index(grace_index_j)
            lmp_index_k = gra_ind.index(grace_index_k)
            lmp_index_l = gra_ind.index(grace_index_l)
            new_txt[0] = str(di_count)
            di_count += 1
            new_txt[2] = str((lmp_index_i + 1) + natm * j)
            new_txt[3] = str((lmp_index_j + 1) + natm * j)
            new_txt[4] = str((lmp_index_k + 1) + natm * j)
            new_txt[5] = str((lmp_index_l + 1) + natm * j)
            new_data.write(" ".join(new_txt))
            new_data.write("\n")

    if (nimp > 0): new_data.write("\nImpropers\n\n")
    imp_count = 1
    for i in range(0, nimp):
        txt = imp_lst[i].split()
        new_tixt = txt[:]
        # if you have more than one molecule, use this loop to create other impropers
        for j in range(0, nmol):
            new_txt = txt[:]
            #put appropriate i and j atom for each impropers
            grace_index_i = int(txt[2]) % (natm + 1) - 1
            grace_index_j = int(txt[3]) % (natm + 1) - 1
            grace_index_k = int(txt[4]) % (natm + 1) - 1
            grace_index_l = int(txt[5]) % (natm + 1) - 1
            lmp_index_i = gra_ind.index(grace_index_i)
            lmp_index_j = gra_ind.index(grace_index_j)
            lmp_index_k = gra_ind.index(grace_index_k)
            lmp_index_l = gra_ind.index(grace_index_l)
            new_txt[0] = str(imp_count)
            imp_count += 1
            new_txt[2] = str((lmp_index_i + 1) + natm * j)
            new_txt[3] = str((lmp_index_j + 1) + natm * j)
            new_txt[4] = str((lmp_index_k + 1) + natm * j)
            new_txt[5] = str((lmp_index_l + 1) + natm * j)
            new_data.write(" ".join(new_txt))
            new_data.write("\n")
    new_data.close()
    os.rename('system.data', 'new_system.data')
    copyfile(os.path.join(grace_dir, 'system.in.init'), 'new_system.in.init')
    copyfile(os.path.join(grace_dir, 'system.in.settings'), 'new_system.in.settings')

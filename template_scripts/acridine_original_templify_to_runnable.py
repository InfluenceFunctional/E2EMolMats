import os

def templify_to_runnable(workdir, original, ori_data_file, New):
    os.chdir(workdir)

    # original = r"4.lt"  # '#input()
    # New = r"5.lt"
    # ori_data_file = r"3.data" #input()

    data = open(original,'r')
    New_data = open(New,'w')
    lines = data.readlines()

    data_file = open(ori_data_file,'r')
    system = open("system.lt",'w')
    system.write("import \"" + New + "\"  # <- defines the \"Acridine\" molecule type.\n\n\n")
    system.write("# Periodic boundary conditions:\nwrite_once(\"Data Boundary\") {\n")
    #print("Is it in orthogonal box?Press 1 if it is. 0 if it is not: ")
    data_file_line = data_file.readline()
    data_file_line = data_file.readline()
    data_file_line = data_file.readline()
    data_file_line = data_file.readline()
    data_file_line = data_file.readline()

    data_file_line = data_file.readline()
    system.write("   " + data_file_line)
    data_file_line = data_file.readline()
    system.write("   " + data_file_line)
    '''xyz box dimensions'''
    data_file_line = data_file.readline()
    system.write("   " + data_file_line)
    data_file_line = data_file.readline()
    system.write("   " + data_file_line)
    data_file_line = data_file.readline()
    system.write("   " + data_file_line)

    data_file_line = data_file.readline()
    if 'xy xz yz' in data_file_line: # if there is a non-orthogonal component, write it
        system.write("   " + data_file_line)

    system.write("}\n\n")
    system.write("# Create 1 \"acridine\" molecules\n# rotated and moved to give polymorph AC1 = new Acridine\n")
    system.write("AC1 = new Acridine")

    counter = 0
    data_atom_checker = 0
    end_data_atom = 0
    data_bond_checker = 0
    end_data_bond = 0
    for ind, i in enumerate(lines):
        counter = counter + 1
        if (counter <5):
            New_data.write(i)
        elif(counter==5):
            New_data.write(i)
            New_data.write("import \"gaff2_acridine.lt\"\nAcridine inherits GAFF2 {\n")
        elif(counter<10):
            New_data.write(i)
        elif(counter < 15):
            line = i.split()
            New_data.write(i)
        elif(data_atom_checker ==0):
            New_data.write(i)
            if(i =="write(\"Data Atoms\") {\n"):
                data_atom_checker =1
        elif(end_data_atom ==0):
            line = i.split()
            if(line[0] =="}"):
                New_data.write(i)
                end_data_atom = 1
            else:
                New_data.write("  " + line[0] +  " " + line[1] + " " + line[2] + " " + line[3] +" " + line[4] + " " + line[5] + " " + line[6] + "\n")
        elif(data_bond_checker ==0):
            if(i =="write(\"Data Bonds\") {\n"):
                data_bond_checker =1
                New_data.write("write(\"Data Bond List\") {\n")
            else:
                New_data.write(i)
        elif(end_data_bond ==0):
            line = i.split()
            if(line[0] =="}"):
                New_data.write(i)
                end_data_bond = 1
            else:
                New_data.write("  " + line[0] + " " + line[2] + " " + line[3] + "\n")
    New_data.write("}")


    New_data.close()
    data.close()
    system.close()



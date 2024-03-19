import os


def moltemp_final(workdir):

    os.chdir(workdir)
    original = "system.data"
    New = "new_" + original
    data = open(original, 'r')
    New_data = open(New, 'w')

    data_file_line = data.readline()  # original atom type. skip
    while data_file_line != "Atoms  # full\n":
        New_data.write(data_file_line)
        data_file_line = data.readline()

    New_data.write("Atoms  # full2\n\n")
#    for i in range(0, 3):
#        data_file_line = data.readline()
#        data_file_line = data.readline()
    # This is right before Masses info is added. I need to change some of ca to ca1 and ca2
    data_file_line = data.readline()
    data_file_line = data.readline()
    mol_id = 1
    atm_counter = 1
    while data_file_line != "\n":
        line = data_file_line.split()
        if atm_counter == 1:
            line[3] = line[3] + " 1"
        elif atm_counter == 2:
            line[3] = line[3] + " 2"
        elif atm_counter == 3:
            line[3] = line[3] + " 3"
        elif atm_counter == 4:
            line[3] = line[3] + " 4"
        elif atm_counter == 5:
            line[3] = line[3] + " 5"
        elif atm_counter == 6:
            line[3] = line[3] + " 6"
        elif atm_counter == 7:
            line[3] = line[3] + " 7"
        elif atm_counter == 8:
            line[3] = line[3] + " 8"
        elif atm_counter == 9:
            line[3] = line[3] + " 9"
        elif atm_counter == 10:
            line[3] = line[3] + " 10"
        elif atm_counter == 11:
            line[3] = line[3] + " 11"
        elif atm_counter == 12:
            line[3] = line[3] + " 12"
        elif atm_counter == 13:
            line[3] = line[3] + " 13"
        elif atm_counter == 14:
            line[3] = line[3] + " 14"
        elif atm_counter == 15:
            line[3] = line[3] + " 15"
        elif atm_counter == 16:
            line[3] = line[3] + " 16"
        elif atm_counter == 17:
            line[3] = line[3] + " 17"
        elif atm_counter == 18:
            line[3] = line[3] + " 18"
        elif atm_counter == 19:
            line[3] = line[3] + " 19"
        elif atm_counter == 20:
            line[3] = line[3] + " 20"
        elif atm_counter == 21:
            line[3] = line[3] + " 21"
        elif atm_counter == 22:
            line[3] = line[3] + " 22"
        elif atm_counter == 23:
            line[3] = line[3] + " 23"
        
        New_data.write(line[0] + " " + line[1] + " " + str(line[2]) + " " + line[3] + " " + line[4] + " " + line[5] + " " + line[6] + "\n");

        # HARD CODED PART FOR BENZAMIDE
        mol_id = line[1]
        data_file_line = data.readline()
        line = data_file_line.split()
        if data_file_line != "\n" and int(mol_id) != int(line[1]):  # in case of banzamide
            atm_counter = 0
        atm_counter += 1

    New_data.write(data_file_line);
    while data_file_line:
        data_file_line = data.readline()
        New_data.write(data_file_line)
    New_data.close()
    data.close()

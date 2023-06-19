import os
os.chdir(r'C:\Users\mikem\crystals\clusters\cluster_structures\test_2')
original = r"1.data"  # '#input()
New = '2.data'  # input()
data = open(original, 'r')
# data = open("nicotinamide_initial_08.data",'r')
New_data = open(New, 'w')
# New_data = open("nicotinamide_python_08.data",'w')
lines = data.readlines()

counter = 0;
mol_counter = 0;
mol = 1;
atm_counter = 1;
for ind,i in enumerate(lines):
    counter = counter + 1;
    if (counter == 3):
        New_data.write("10 atom types\n")
    elif (counter < 9):
        New_data.write(i);
    elif (counter == 11):
        New_data.write("\nMasses\n\n1 1.008  # ha\n2 1.008  # h4\n3 1.008  # hn\n4 14.01  # n\n5 14.01  # nb\n6 12.01  # c\n7 12.01  # ca\n8 12.01  # ca1\n9 12.01  # ca2\n10 16.0  # o\n")
        New_data.write(i);
    elif (counter == 17):
        New_data.write(i)
        New_data.write("\n")
    elif counter>=19:
        line = i.split();
        if (mol_counter == 15):
            mol = mol + 1;
            mol_counter = 0;
        mol_counter = mol_counter + 1
        line[1] = mol
        # print(line[1])
        if atm_counter == 1:
            line[2] = 8
            line[3] = -0.326615;
            atm_counter = atm_counter + 1;
        elif atm_counter == 2:
            line[2] = 7;
            line[3] = 0.380568;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 3):
            line[2] = 5;
            line[3] = -0.585364;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 4):
            line[2] = 9;
            line[3] = 0.384166;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 5):
            line[2] = 7;
            line[3] = -0.38538;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 6):
            line[2] = 7;
            line[3] = 0.173788;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 7):
            line[2] = 6;
            line[3] = 0.6054;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 8):
            line[2] = 10;
            line[3] = -0.479634;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 9):
            line[2] = 4;
            line[3] = -0.779885;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 10):
            line[2] = 3;
            line[3] = 0.357505;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 11):
            line[2] = 3;
            line[3] = 0.357505;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 12):
            line[2] = 2;
            line[3] = 0.027993;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 13):
            line[2] = 2;
            line[3] = 0.034858;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 14):
            line[2] = 1;
            line[3] = 0.157212;
            atm_counter = atm_counter + 1;
        elif (atm_counter == 15):
            line[2] = 1;
            line[3] = 0.077882;
            atm_counter = atm_counter + 1;
        if (atm_counter == 16):
            atm_counter = 1;
        for j in range(0, 7):
            New_data.write(str(line[j]) + ' ')
        New_data.write("\n")

New_data.close()
data.close()

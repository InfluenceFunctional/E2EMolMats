print("Enter the file name you just created: ")
original = input()
New = "new_" + original
data = open(original,'r')
#data = open("nicotinamide_initial_08.data",'r')
New_data = open(New,'w')
#New_data = open("nicotinamide_python_08.data",'w')
lines = data.readlines()

print("Enter the the name of the data file: ")
ori_data_file = input()
data_file = open(ori_data_file,'r');
system = open("system.lt",'w');
system.write("import \"" + New + "\"  # <- defines the \"Urea\" molecule type.\n\n\n")
system.write("# Periodic boundary conditions:\nwrite_once(\"Data Boundary\") {\n")
ortho = 0;
print("Is it in orthogonal box?Press 1 if it is. 0 if it is not: ")
ortho = input();
data_file_line = data_file.readline()
data_file_line = data_file.readline()
data_file_line = data_file.readline()
data_file_line = data_file.readline()
data_file_line = data_file.readline()

data_file_line = data_file.readline()
system.write("   " + data_file_line)
data_file_line = data_file.readline()
system.write("   " + data_file_line)
data_file_line = data_file.readline()
system.write("   " + data_file_line)
if(ortho =="1"):
    data_file_line = data_file.readline()
    system.write("   " + data_file_line)
system.write("}\n\n")
system.write("# Create 1 \"nicotinamide\" molecules\n# rotated and moved to give polymorph Inic1 = new Nicotinamide\n")
system.write("nic1 = new Nicotinamide")

counter = 0;
data_atom_checker = 0;
end_data_atom = 0;
data_bond_checker = 0;
end_data_bond = 0;
for i in lines:
    counter = counter + 1;
    if (counter <5):
        New_data.write(i);
    elif(counter==5):
        New_data.write(i);
        New_data.write("import \"gaff2_nicotinamid_nolong.lt\"\nNicotinamide inherits GAFF2 {\n");
    elif(counter<10):
        New_data.write(i);
    elif(counter < 22):
        line = i.split();
        if(line[0] =="@atom:ca1" or line[0] =="@atom:ca2"):
            continue;
        else:
            New_data.write(i);
    elif(data_atom_checker ==0):
        New_data.write(i);
        if(i =="write(\"Data Atoms\") {\n"):
            data_atom_checker =1;
    elif(end_data_atom ==0):
        line = i.split();
        if(line[0] =="}"):
            New_data.write(i);
            end_data_atom = 1;
        else:
            if(line[2] =="@atom:ca1" or line[2] =="@atom:ca2"):
                line[2] = "@atom:ca";
            New_data.write("  " + line[0] +  " " + line[1] + " " + line[2] + " " + line[3] +" " + line[4] + " " + line[5] + " " + line[6] + "\n")
    elif(data_bond_checker ==0):
        if(i =="write(\"Data Bonds\") {\n"):
            data_bond_checker =1;
            New_data.write("write(\"Data Bond List\") {\n")
        else:
            New_data.write(i);
    elif(end_data_bond ==0):
        line = i.split();
        if(line[0] =="}"):
            New_data.write(i);
            end_data_bond = 1;
        else:
            New_data.write("  " + line[0] + " " + line[2] + " " + line[3] + "\n")
New_data.write("}")


New_data.close()
data.close()
system.close()

    

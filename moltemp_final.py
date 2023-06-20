original = "system.data"
New = "new_"+original
data = open(original,'r')
New_data = open(New,'w')

for i in range(0,8):
    data_file_line = data.readline()
    New_data.write(data_file_line)
data_file_line = data.readline()#original atom type. skip
New_data.write("     10  atom types\n")
while(data_file_line !="Masses\n"):
    data_file_line = data.readline()
    New_data.write(data_file_line)

for i in range(0,9):
    data_file_line = data.readline()
    New_data.write(data_file_line)
New_data.write("9 12.01 # ca1\n10 12.01 # ca2\n")
for i in range(0,3):
    data_file_line = data.readline()
    New_data.write(data_file_line)
#This is right before Masses info is added. I need to change some of ca to ca1 and ca2
data_file_line = data.readline()
while(data_file_line !="\n"):
    line = data_file_line.split()
    if(int(line[0])%15==1):
        line[2]=9;
    elif(int(line[0])%15==4):
        line[2]=10;
    New_data.write(line[0] + " " + line[1] + " " + str(line[2]) + " " + line[3] + " " + line[4] + " " + line[5] + " " + line[6] + "\n");
    data_file_line = data.readline();
New_data.write(data_file_line);
while(data_file_line):
    data_file_line = data.readline()
    New_data.write(data_file_line)
New_data.close();
data.close();

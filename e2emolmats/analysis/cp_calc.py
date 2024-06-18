import os
import math

''' MK overall notes
- semicolons are no longer necessary in python. --DK: Yes, I do know that, but I just have habit of writing one since I do write in C++ as well which require semicolons.
- not using PEP8 spacing rules. Not a big deal just a note. e.g., --DK: Thanks for pointing out this. I think I was not consistant about spacing.
    # ['abc', '123']
    instead of
    #['abc','123']
- variable names are not always clear   --DK: Sorry about this. I tried not to change the name aster I assign to reduce the potential error.
- calculations need to be referenced or at least commented.     --DK: Yes, I sent slack maeeage with equation and reference textbook.
- calculations should also be broken into simpler steps, currently it is very hard to follow, why? We can make it very easy if we want. --DK: Yes, I think I can call more variable to substitute the complicated calculation.
'''

dir_path = os.path.dirname(os.path.realpath(__file__))
#seedLst = ["seed1","seed2","seed3","seed4","seed5","seed6","seed7","seed8","seed9","seed10",]
seedLst = ["seed3", "seed4", "seed5", "seed6", "seed7"]
#seedLst = ["seed10"]
#seedLst = ["seed3"]
TempLst = ["T100", "T300"]
#TempLst = ["T10","T100","T300","T400","T450","T500","T550","T600"]
#TempLst = ["T100"]
#TempLst = ["T300","T400","T450","T500","T550","T600"]
crystal = [2, 3]
#crystal = [2]
#nmol = [336,504,504,448,448,384,384]
nmol = [336, 504]
natom = 7728

# MK confirm these and comment the units
kb = 0.008314462618  #boltzman factor. checked(https://en.wikipedia.org/wiki/Boltzmann_constant). kJ/mol/K
kJ = 4.184  #kJ/kcal
nktv2p = 68568.415  # MK what is this?  #convert atm*Ang^3 to kcal/mol. From LAMMPS. (kcal/mol)/(atm*Ang^3)?
aveElst = []
aveE2lst = []
aveentlst = []
aveent2lst = []
avePElst = []
avePE2lst = []
aveKElst = []
aveKE2lst = []
aveTmplst = []
aveCvlst = []
aveCpElst = []
aveCvlst_KE = []
P = 1
for seed in seedLst:
    for temperature in TempLst:
        # MK what is this? are we recording the temperature? I think you can do int(string) which might be simpler
        #DK: reading TempLst I guess I don't need this anymore.
        print(temperature)
        for j in crystal:
            directory = "form" + str(j)
            init = open(seed + "/" + temperature + "/" + directory + "/log.lammps", 'r')
            after = open(seed + "/" + temperature + "/" + directory + "/log_after.lammps", 'w')
            #            run_ave = open(l + "/" + k + "/" + directory + "/log_run_ave.lammps", 'w')
            #            run_ave.write("   Step          Time           <H>2           <H2>               Cp             Cperror\n")
            test = open(seed + "/" + temperature + "/" + directory + "/test.lammps", 'w')

            # MK surely there's a more efficient way to do this initialization?
            # maybe with a dictcomp?
            #DK: Sure I think we can try that too.
            before_log = True
            after_log = False
            counter = 0
            tmp = 0
            entlst = []
            ent2lst = []
            Elst = []
            E2lst = []
            PElst = []
            PE2lst = []
            KElst = []
            KE2lst = []
            entlst = []
            ent2lst = []
            tmplst = []
            ent_ave = 0
            ent2_ave = 0
            E_ave = 0
            E2_ave = 0
            PE_ave = 0
            PE2_ave = 0
            KE_ave = 0
            KE2_ave = 0
            ent_ave = 0
            ent2_ave = 0
            tmp_ave = 0
            err1 = 0
            err2 = 0
            errE1 = 0
            errE2 = 0
            final = 0
            mol = nmol[crystal.index(j)]
            for i in init:
                if (before_log == False and after_log == False):
                    #print(i)
                    line = i.split()
                    after.write(i)
                    counter += 1
                    if (len(line) != 17):
                        after_log = True
                        counter -= 1
                        continue
                    # MK it is very hard to follow the calculations at the same time as you are parsing the text file
                    # I highly recomment extracting all the relevant information first, and then doing calculations after
                    # or at minimum, assigning all these float(line[ind]) values to variables so we can confirm what we are
                    # looking at

                    #DK: Thanks for the suggestion. Maybe I can extract information in this loop and calculate after the loop.

                    tmp += float(line[2])

                    #enthalpy = total energy + PV where P is constant and V is the instantanious value.
                    #P *  float(line[4]) / nktv2p means PV in kcal/mol
                    Elst.append((float(line[5]) + P * float(line[4]) / nktv2p) * kJ)
                    E2lst.append((((float(line[5]) + P * float(line[4]) / nktv2p) * kJ) ** 2))
                    entlst.append(float(line[6]) * kJ)
                    ent2lst.append(((float(line[6]) * kJ) ** 2))
                    tmplst.append(float(line[2]))
                    test.write((line[6]) + " " + line[2] + " " + line[3] + " " + line[4] + " \n")

                if ("Loop time" in i):
                    after_log = True
                if (("   Step          Time") in i):
                    after.write(i)
                    before_log = False
                    after_log = False
            init.close()
            after.close()
            #run_ave.close()
            #ent_ave = ent/counter
            #ent2_ave = ent2/counter

            first_half = int(len(entlst) / 2)
            second_half = len(entlst) - first_half
            #I tried to get the average of energy. I take the average of the second half of the simulation (for example if I have 10ns simulation, I take the average from 5-10ns)
            for m in range(second_half, len(entlst)):  # MK no idea what's going on here
                ent_ave += entlst[m] / second_half
                ent2_ave += ent2lst[m] / second_half
                E_ave += Elst[m] / second_half
                E2_ave += E2lst[m] / second_half
                tmp_ave += tmplst[m] / second_half

            tmp_ave = tmp / counter
            test.close()
            aveentlst.append(ent_ave)
            aveent2lst.append(ent2_ave)
            aveElst.append(E_ave)
            aveE2lst.append(E2_ave)
            aveTmplst.append(tmp_ave)
            aveCvlst.append(0)
            aveCpElst.append(0)
            aveCvlst_KE.append(0)

            init.close()
            after.close()
        print(counter)

#average of average. vvvvv From here, I tried to calculate the average over seeds. Tried to get the "final" result.  vvvvv
ave_of_ave_tmp = []
ave_of_ave_ent = []
ave_of_ave_ent2 = []
ave_of_ave_E = []
ave_of_ave_E2 = []
ave_of_ave_PE = []
ave_of_ave_PE2 = []
ave_of_ave_KE = []
ave_of_ave_KE2 = []
err1 = []
err2 = []
errE1 = []
errE2 = []
err1KE = []
err2KE = []
cov_each = []
covE_each = []
cov_eachKE = []

sigmaent = []
sigmaent2 = []
sigmaE = []
sigmaE2 = []
sigmaPE = []
sigmaPE2 = []
sigmaKE = []
sigmaKE2 = []
cov_ave = []
covE_ave = []
cov_aveKE = []
stdCv = []
stdCpE = []
stdCv_KE = []

Cv = []
CpE = []
Cv_KE = []
seed_num = len(seedLst)
temp_num = len(TempLst)
crystal_num = len(crystal)

#assigning from each temp and polymorph
for j in range(0, temp_num * crystal_num):
    ave_of_ave_tmp.append(0)
    ave_of_ave_ent.append(0)
    ave_of_ave_ent2.append(0)
    ave_of_ave_E.append(0)
    ave_of_ave_E2.append(0)
    ave_of_ave_PE.append(0)
    ave_of_ave_PE2.append(0)
    ave_of_ave_KE.append(0)
    ave_of_ave_KE2.append(0)
    Cv.append(0)
    CpE.append(0)
    Cv_KE.append(0)
    err1.append(0)
    err2.append(0)
    errE1.append(0)
    errE2.append(0)
    err1KE.append(0)
    err2KE.append(0)
    cov_each.append(0)
    covE_each.append(0)
    cov_eachKE.append(0)
    sigmaent.append(0)
    sigmaent2.append(0)
    sigmaE.append(0)
    sigmaE2.append(0)
    sigmaPE.append(0)
    sigmaPE2.append(0)
    sigmaKE.append(0)
    sigmaKE2.append(0)
    cov_ave.append(0)
    covE_ave.append(0)
    cov_aveKE.append(0)
    stdCv.append(0)
    stdCpE.append(0)
    stdCv_KE.append(0)

#Loop over all the information of temprature, ent, and E and get the average over seeds
for i in range(0, seed_num):
    for temperature in range(0, temp_num):
        for j in range(0, crystal_num):
            #each_crystal_temp is the index of each crystal and temperature
            #whole_loop is each_crystal_temp with coresponding seed
            each_crystal_temp = temperature * crystal_num + j
            whole_loop = i * temp_num * crystal_num + temperature * crystal_num + j

            ave_of_ave_tmp[each_crystal_temp] += aveTmplst[whole_loop] / seed_num
            ave_of_ave_ent[each_crystal_temp] += aveentlst[whole_loop] / seed_num
            ave_of_ave_ent2[each_crystal_temp] += aveent2lst[whole_loop] / seed_num
            ave_of_ave_E[each_crystal_temp] += aveElst[whole_loop] / seed_num
            ave_of_ave_E2[each_crystal_temp] += aveE2lst[whole_loop] / seed_num

#Here we calculate error and covariance using the average ent and E calculated above
for i in range(0, seed_num):
    for temperature in range(0, temp_num):
        for j in range(0, crystal_num):
            #each_crystal_temp is the index of each crystal and temperature
            #whole_loop is each_crystal_temp with coresponding seed
            each_crystal_temp = temperature * crystal_num + j
            whole_loop = i * temp_num * crystal_num + temperature * crystal_num + j

            err1[each_crystal_temp] += ((aveentlst[whole_loop] - ave_of_ave_ent[each_crystal_temp]) ** 2)
            err2[each_crystal_temp] += ((aveent2lst[whole_loop] - ave_of_ave_ent2[each_crystal_temp]) ** 2)
            errE1[each_crystal_temp] += ((aveElst[whole_loop] - ave_of_ave_E[each_crystal_temp]) ** 2)
            errE2[each_crystal_temp] += ((aveE2lst[whole_loop] - ave_of_ave_E2[each_crystal_temp]) ** 2)
            cov_each[each_crystal_temp] += (aveentlst[whole_loop] - ave_of_ave_ent[each_crystal_temp]) * (
                    aveent2lst[whole_loop] - ave_of_ave_ent2[each_crystal_temp])
            covE_each[each_crystal_temp] += (aveElst[whole_loop] - ave_of_ave_E[each_crystal_temp]) * (
                    aveE2lst[whole_loop] - ave_of_ave_E2[each_crystal_temp])

#calculate the std
for j in range(0, temp_num * crystal_num):
    sigmaent[j] = math.sqrt(err1[j] / (seed_num))
    sigmaent2[j] = math.sqrt(err2[j] / (seed_num))
    sigmaE[j] = math.sqrt(errE1[j] / (seed_num))
    sigmaE2[j] = math.sqrt(errE2[j] / (seed_num))
    cov_ave[j] = cov_each[j] / (seed_num)
    covE_ave[j] = covE_each[j] / (seed_num)

#Calculate Cp or Cv
for temperature in range(0, temp_num):
    for j in range(0, crystal_num):
        # MK with serious math here, we definitely need some reference to tell us what is going on, or at least, we should break it into steps which we can debug line-by-line.
        # Currently, each line is way too complex to debug at a glance. We should try to make it very simple and easy
        # I agree. I did the test calculation before in different python code. We can still try with having only one seed, one temperature, for one polymorph. I have google seet that calculate the sample.
        #(This is for calculating Cv. https://docs.google.com/spreadsheets/d/1wLfvugIS9j6RsKS6KoqZB65rTDlzOlB5kg67YU-32Zs/edit?usp=sharing)

        #each_crystal_temp is the index of each crystal and temperature
        each_crystal_temp = temperature * crystal_num + j
        kbTT = kb * ave_of_ave_tmp[each_crystal_temp] * ave_of_ave_tmp[each_crystal_temp]

        #RMS = <dE^2> = <E^2> - <E>^2
        RMS_ent = (-1.0 * ave_of_ave_ent[each_crystal_temp] * ave_of_ave_ent[each_crystal_temp] + ave_of_ave_ent2[
            each_crystal_temp])
        RMS_E = (-1.0 * ave_of_ave_E[each_crystal_temp] * ave_of_ave_E[each_crystal_temp] + ave_of_ave_E2[
            each_crystal_temp])

        Cv[each_crystal_temp] = RMS_ent / kbTT / nmol[j]
        CpE[each_crystal_temp] = RMS_E / kbTT / nmol[j]

        dCv_dent2_sigma2 = sigmaent2[each_crystal_temp] / kbTT
        dCv_dent_sigma = -2 * ave_of_ave_ent[each_crystal_temp] * sigmaent[each_crystal_temp] / kbTT
        cov_ent_ent2 = 2 * (-2 * ave_of_ave_ent[each_crystal_temp] * cov_ave[each_crystal_temp] / kbTT / kbTT)

        stdCv[each_crystal_temp] = math.sqrt((dCv_dent2_sigma2) ** 2 + (dCv_dent_sigma) ** 2 + cov_ent_ent2) / nmol[j]

        dCpE_dE2_sigma2 = sigmaE2[each_crystal_temp] / kbTT
        dCpE_dE_sigma = -2 * ave_of_ave_E[each_crystal_temp] * sigmaE[each_crystal_temp] / kbTT
        cov_E_E2 = 2 * (-2 * ave_of_ave_E[each_crystal_temp] * covE_ave[each_crystal_temp] / kbTT / kbTT)

        stdCpE[each_crystal_temp] = math.sqrt((dCpE_dE2_sigma2) ** 2 + (dCpE_dE_sigma) ** 2 + cov_E_E2) / nmol[j]

#calculate Cv and Cp. no average.
for i in range(0, seed_num):
    for temperature in range(0, temp_num):
        for j in range(0, crystal_num):
            #each_crystal_temp is the index of each crystal and temperature
            #whole_loop is each_crystal_temp with coresponding seed
            each_crystal_temp = temperature * crystal_num + j
            whole_loop = i * temp_num * crystal_num + temperature * crystal_num + j

            A = (-1.0 * aveentlst[whole_loop] * aveentlst[whole_loop] + aveent2lst[whole_loop])
            A_E = (-1.0 * aveElst[whole_loop] * aveElst[whole_loop] + aveE2lst[whole_loop])

            aveCvlst[whole_loop] = A / kb / aveTmplst[whole_loop] / aveTmplst[whole_loop] / nmol[j]
            aveCpElst[whole_loop] = A_E / kb / aveTmplst[whole_loop] / aveTmplst[whole_loop] / nmol[j]

print("##################################################################")
print("ave_tmp          aveCpElst")
for i in range(0, crystal_num * temp_num * seed_num):
    print(str(aveTmplst[i]) + " " + str(aveCpElst[i]))
print("##################################################################")

print("crystal ave_tmp      CpE                          stdCpE")
print(len(CpE))
print(crystal_num * temp_num)
for j in range(0, crystal_num * temp_num):
    print(str(crystal[j % crystal_num]) + " " + str(ave_of_ave_tmp[j]) + " " + str(CpE[j]) + " " + str(stdCpE[j]))
print(ave_of_ave_E)
print(ave_of_ave_E2)
print(CpE)

import os
from distutils.dir_util import copy_tree
from shutil import copyfile
from time import sleep
from random import shuffle

classifier_source_path = '/scratch/mk8347/daisuke_classifier/molvec_lib/'
run_path = '/scratch/mk8347/molecule_clusters/defect_clusters_5/'
os.chdir(run_path)

run_dirs = os.listdir()
#
# os.system('module load intel/19.1.2')
# os.system('module load openmpi/intel/4.1.1')

shuffle(run_dirs)  # make parallel runs feasible

for dir_path in run_dirs[:5]:
    try_iter = False
    finished = False
    wait_steps = 0
    try:
        print(int(dir_path))
        try_iter = True
    except:
        pass

    if try_iter and (not os.path.exists(dir_path + '/run_ave_traj.dump')):
        copy_tree(classifier_source_path, dir_path)
        copyfile('example/weight_nicotinamide.txt', dir_path + 'weight.txt')
        copyfile('example/INPUT_nico_cluster2.dat', dir_path + 'INPUT.dat')
        os.chdir(dir_path)
        os.system('make')
        sleep(20)
        os.system(f'./sa traj.dump lammps_vec INPUT.dat')

        while (not finished) or (wait_steps < 20):  # wait a maximum of 200 seconds before giving up
            if os.path.exists('run_ave_traj.dump'):
                finished = True
            else:
                sleep(10)
                wait_steps += 1

        os.chdir('../../')


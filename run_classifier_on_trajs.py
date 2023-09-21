import os
from distutils.dir_util import copy_tree
from shutil import copyfile
from time import sleep

classifier_source_path = '/scratch/mk8347/daisuke_classifier/molvec_lib/'
run_path = '/scratch/mk8347/molecule_clusters/dev10/'
os.chdir(run_path)

run_dirs = os.listdir()

os.system('module load intel/19.1.2')
os.system('module load openmpi/intel/4.1.1')

copy_tree(classifier_source_path, './')
copyfile('example/weight_nicotinamide.txt','weight.txt')
copyfile('example/INPUT_nico_cluster2.dat', 'INPUT.dat')

os.system('make')

for dir_path in run_dirs:
    try_iter = False
    finished = False
    wait_steps = 0
    try:
        print(int(dir_path))
        try_iter = True
    except:
        pass

    if try_iter:
        os.system(f'./sa {dir_path}/traj.dump lammps_vec INPUT.dat')
        sleep(10)

        while (not finished) or (wait_steps < 20): # wait a maximum of 200 seconds before giving up
            if os.path.exists('run_ave_traj.dump'):
                finished = True
                os.rename('run_ave_traj.dump', dir_path + '/run_ave_traj.dump')
                os.rename('new_traj.dump', dir_path + '/new_traj.dump')
            else:
                sleep(10)
                wait_steps += 1


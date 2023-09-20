import os
from distutils.dir_util import copy_tree
from shutil import copyfile

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
    try:
        print(int(dir_path))
        copyfile(dir_path + 'traj.dump', './')
        os.system('.sa traj.dump lammps_vec INPUT.dat')
        os.system('rm traj.dump')
        os.rename('run_ave_traj.dump', dir_path + 'run_ave_traj.dump')
        os.rename('new_traj.dump', dir_path + 'new_traj.dump')
    except:
        pass

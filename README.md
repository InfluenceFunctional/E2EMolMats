convert xyz files of molecule clusters to LAMMPS inputs
nicotinamide only but with options for vacancies, random reorientations and benzamide defects (last will not work yet for all purposes)

edit paths in generate_cluster_structures and xyz_to_lammps to your local specs, then run xyz_to_lammps

requires ase, numpy, moltemplate, ovito, scipy

"""run configs to loop over - will take all combinations of the below lists - grows combinatorially!!"""
"""run this script to make the corresponding .yaml file, which is used at runtime"""
import yaml

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[10, 10, 10]],  # size of initial bulk supercell, from which finite subsamples may be carved.
    'temperature': [100],  # Kelvin - the temperature at which to carry out primary sampling
    'structure_identifier': ['acridine/Form2'],
    # ['acridine/Form2'], # ['nicotinamide/NICOAM17'], #'["nicotinamide/NICOAM16"],  # , "NICOAM17"],  # names of molecule and polymorph structures to be used as bases.
    'defect_rate': [.5],
    # fraction of molecules to be substituted with appropriately aligned defects - only works for Benzamide in Nicotinamide
    'defect_type': ['2,7-dihydroxynaphthalene'],
    # what molecule to substitute in the lattice - 'benzamide' or 'isonicotinamide' for nicotinamide, 'anthracene' or '2,7-dihydroxynaphthalene' for acridine
    'gap_rate': [0],  # fraction of molecule sites to be left vacant
    'scramble_rate': [0],  # fraction of molecules to be randomly rotated
    'seed': [1],  # integers greater than zero
    'damping': [str(100.0)],  # for LAMMPS Langevin dynamics
    'max_sphere_radius': [20],  # if carving a finite cluster from a bulk structure, the radius of the sphere

    # static items - DO NOT SET AS LIST
    'cluster_type': 'supercell',
    # type of structure to simulate. "supercell" a nxnxn bulk crystal supercell. "spherical" a finite cluster in vacuum.
    'run_time': int(1e5),  # sampling time in femtoseconds
    'box_type': 'p',
    # box type in LAMMPS dimensions 'p' for periodic typically used even for vacuum simulations, just with very large box
    'integrator': 'nosehoover',  # nosehoover, npt, nvt
    'pressure_direction': 'z',  # iso, x, y, or z for npt pressure directionality. Also the interface axis in interfaces.
    'ramp_temperature': False,  # linearly ramp temperature in main sampling run from 0-temperature
    'init_temperature': 0,
    'print_steps': int(1e2),  # how many timepoints to print in sampling trajectory
    'min_inter_cluster_distance': 5,  # 40,  # sets periodic box size in cluster simulations
    'bulk_crystal': True,  # if true, periodic structure
    'machine': 'local',  # 'local' or 'cluster' have different associated paths
    'run_name': 'dev',
    'min_lattice_length': 5,
    # for periodic bulk simulations. Supercell a x b x c a,b,c will be set to approximately at least this edge length.
    'prep_crystal_in_melt': True,  # prepare a frozen nanocrystal in a melted environment
    'prep_melt_interface': False,  # split supercell in half and melt half then reequilibrate
    'prep_bulk_melt': False,  # prepare a bulk melted structure - npt equil, nvt melt, nvt cool, npt equil
    'equil_time': 1e5,  # equilibration time, mostly for prep_crystal_in_melt steps
    'melt_temperature': 2000,  # melt temperature of prep_crystal_in_melt runs
    'atom_style': 'full',  # 'full' or 'full2' depending on if we want symmetry information
    'submit_lammps_slurm': False,
}

filename = __file__
with open(filename.split('.py')[0] + '.yaml', 'w') as outfile:
    yaml.dump(batch_config, outfile, default_flow_style=False)

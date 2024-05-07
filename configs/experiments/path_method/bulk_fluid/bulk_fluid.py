"""run configs to loop over - will take all combinations of the below lists - grows combinatorially!!"""
"""run this script to make the corresponding .yaml file, which is used at runtime"""
import yaml

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[15, 15, 15]],  # size of initial bulk supercell, from which finite subsamples may be carved.
    'temperature': [700],  # Kelvin - the temperature at which to carry out primary sampling
    'structure_identifier': ['acridine/Form4'],
    # ['acridine/Form2'], # ['nicotinamide/NICOAM17'], #'["nicotinamide/NICOAM16"],  # , "NICOAM17"],  # names of molecule and polymorph structures to be used as bases.
    'defect_rate': [0],
    # fraction of molecules to be substituted with appropriately aligned defects - only works for Benzamide in Nicotinamide
    'defect_type': ['anthracene'],
    # what molecule to substitute in the lattice - 'benzamide' or 'isonicotinamide' for nicotinamide, 'anthracene' or '2,7-dihydroxynaphthalene' for acridine
    'gap_rate': [0],  # fraction of molecule sites to be left vacant
    'scramble_rate': [0],  # fraction of molecules to be randomly rotated
    'seed': [1],  # integers greater than zero
    'damping': [str(100.0)],  # for LAMMPS Langevin dynamics
    'max_sphere_radius': [10],  # if carving a finite cluster from a bulk structure, the radius of the sphere

    # static items - DO NOT SET AS LIST
    'cluster_type': 'supercell',
    # type of structure to simulate. "supercell" a nxnxn bulk crystal supercell. "spherical" a finite cluster in vacuum.
    'run_time': 2000000,  # sampling time in femtoseconds
    'box_type': 'p',
    # box type in LAMMPS dimensions 'p' for periodic typically used even for vacuum simulations, just with very large box
    'integrator': 'npt_iso',  # nosehoover, npt_iso, npt_aniso, npt_tri, nvt
    'fix_com': False,  # fix center of mass motion
    'ramp_temperature': True,  # linearly ramp temperature in main sampling run from 0-temperature
    'init_temperature': 700,
    'print_steps': int(2e2),  # how many timepoints to print in sampling trajectory
    'min_inter_cluster_distance': 0,  # 40,  # sets periodic box size in cluster simulations
    'bulk_crystal': True,  # if true, periodic structure
    'machine': 'cluster',  # 'local' or 'cluster' have different associated paths
    'run_name': 'bulk_fluid',
    'min_lattice_length': 40,
    # for periodic bulk simulations. Supercell a x b x c a,b,c will be set to approximately at least this edge length.
    'prep_crystal_in_melt': False,  # Work in progress - prepare a frozen nanocrystal in a melted environment
    'prep_melt_interface': False,  # Work in progress - split supercell in half along the fractional z direction
    'prep_bulk_melt': False,  # prepare a bulk melted structure - npt equil, nvt melt, nvt cool, npt equil
    'equil_time': 2000000,  # equilibration time, mostly for prep_crystal_in_melt steps
    'melt_temperature': 2000,  # melt temperature of prep_crystal_in_melt runs
    'atom_style': 'full',  # 'full' or 'full2' depending on if we want symmetry information
    'submit_lammps_slurm': False,
}

filename = __file__
with open(filename.split('.py')[0] + '.yaml', 'w') as outfile:
    yaml.dump(batch_config, outfile, default_flow_style=False)

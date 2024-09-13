"""run configs to loop over - will take all combinations of the below lists - grows combinatorially!!"""
"""run this script to make the corresponding .yaml file, which is used at runtime"""
"""this test is to see if the new force field is working / general functionality"""
import yaml

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[10, 10, 10]],  # size of initial bulk supercell, from which finite subsamples may be carved. Should generally be *very large*
    'temperature': [330, 340, 350, 360, 370, 380, 390, 400, 410],  # Kelvin
    'structure_identifier': ['acridine/Form2',
                             'acridine/Form4'],
    'defect_rate': [0],
    'defect_type': [None],
    # what molecule to substitute in the lattice - 'benzamide' or 'isonicotinamide' for nicotinamide, 'anthracene' or '2,7-dihydroxynaphthalene' for acridine
    'gap_rate': [0],  # fraction of molecule sites to be left vacant
    'scramble_rate': [0],  # fraction of molecules to be randomly rotated
    'seed': [1, 2],  # integers greater than zero
    'damping': [str(100.0)],  # for LAMMPS Langevin dynamics
    'max_sphere_radius': [10],  # if carving a finite cluster from a bulk structure, the radius of the sphere
    'invert_defects': False,

    # static items - DO NOT SET AS LIST
    'run_time': 1e8,  # sampling time in femtoseconds
    'cluster_type': 'supercell',
    # type of structure to simulate. "supercell" a nxnxn bulk crystal supercell. "spherical" a finite cluster in vacuum.
    'box_type': 'p',
    # box type in LAMMPS dimensions 'p' for periodic typically used even for vacuum simulations, just with very large box
    'integrator': 'npt',  # nosehoover, npt, nvt
    'pressure_direction': ['x', 'y', 'z'],  # iso, x, y, or z for npt pressure directionality
    'ramp_temperature': False,  # linearly ramp temperature in main sampling run from 0-temperature
    'init_temperature': 200,  # for ramps only
    'print_steps': int(1e2),  # how many timepoints to print in sampling trajectory
    'min_inter_cluster_distance': 0,  # 40,  # sets periodic box size in cluster simulations, 0 or None if unused
    'bulk_crystal': True,  # if true, periodic structure
    'machine': 'cluster',  # 'local' or 'cluster' have different associated paths
    'run_name': 'acridine_melt_interface2',
    'min_lattice_length': [25],
    # for periodic bulk simulations. Supercell a x b x c a,b,c will be set to approximately at least this edge length.
    'prep_crystal_in_melt': False,  # Work in progress - prepare a frozen nanocrystal in a melted environment
    'prep_melt_interface': True,  # Work in progress - split supercell in half along the fractional z direction
    'prep_bulk_melt': False,  # prepare a bulk melted structure - npt equil, nvt melt, nvt cool, npt equil
    'equil_time': 5e5,  # equilibration time, for melt preparation steps
    'melt_temperature': 2000,  # melt temperature of prep_crystal_in_melt runs
    'atom_style': 'full',  # 'full' or 'full2' depending on if we want symmetry information
    'submit_lammps_slurm': False,
}

filename = __file__
with open(filename.split('.py')[0] + '.yaml', 'w') as outfile:
    yaml.dump(batch_config, outfile, default_flow_style=False)

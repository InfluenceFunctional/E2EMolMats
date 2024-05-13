"""run configs to loop over - will take all combinations of the below lists - grows combinatorially!!"""

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[15, 15, 15]],  # size of initial bulk supercell, from which finite subsamples may be carved.
    'temperature': [350],  # Kelvin
    'structure_identifier': ['nicotinamide/NICOAM17'],  # Form IX
    'defect_rate': [0],  # fraction of molecules to be substituted with appropriately aligned defects - only works for Benzamide in Nicotinamide
    'defect_type': ['benzamide'],  # what molecule to substitute in the lattice - 'benzamide' or 'isonicotinamide' for nicotinamide, 'anthracene' or '2,7-dihydroxynaphthalene' for acridine
    'gap_rate': [0],  # fraction of molecule sites to be left vacant
    'scramble_rate': [0],  # fraction of molecules to be randomly rotated
    'seed': [1, 2],
    'damping': [str(100.0)],  # for LAMMPS Langevin dynamics
    'max_sphere_radius': [10, 12.5, 15, 17.5, 20, 22.5, 25],  # if carving a finite cluster from a bulk structure, the radius of the sphere

    # static items - DO NOT SET AS LIST
    'cluster_type': 'supercell',
    'run_time': int(1e8),
    'box_type': 'p',
    'integrator': 'nosehoover',
    'print_steps': int(1e3),
    'min_inter_cluster_distance': 20,
    'bulk_crystal': True,
    'machine': 'cluster',
    'run_name': 'crystal_in_melt_test11',
    'min_lattice_length': 10,
    'prep_crystal_in_melt': True,
    'equil_time': int(2.5e5),
    'melt_temperature': 2000,
}

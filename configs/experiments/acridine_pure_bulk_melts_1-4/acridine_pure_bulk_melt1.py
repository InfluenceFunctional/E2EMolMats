"""run configs to loop over - will take all combinations of the below lists - grows combinatorially!!"""

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[15, 15, 15]],  # size of initial bulk supercell, from which finite subsamples may be carved.
    'temperature': [300, 350, 400],  # Kelvin
    'structure_identifier': ['acridine/Form4'],  #
    'defect_rate': [0],  # fraction of molecules to be substituted with appropriately aligned defects - only works for Benzamide in Nicotinamide
    'defect_type': ['anthracene'],  # what molecule to substitute in the lattice - 'benzamide' or 'isonicotinamide' for nicotinamide, 'anthracene' or '2,7-dihydroxynaphthalene' for acridine
    'gap_rate': [0, 0.1, 0.2],  # fraction of molecule sites to be left vacant
    'scramble_rate': [0],  # fraction of molecules to be randomly rotated
    'seed': [1],
    'damping': [str(100.0)],  # for LAMMPS Langevin dynamics
    'max_sphere_radius': [None],  # if carving a finite cluster from a bulk structure, the radius of the sphere

    # static items - DO NOT SET AS LIST
    'cluster_type': 'supercell',
    'run_time': int(1e7),
    'box_type': 'p',
    'integrator': 'npt',
    'print_steps': int(1e3),
    'min_inter_cluster_distance': 0,
    'bulk_crystal': True,
    'machine': 'cluster',
    'run_name': 'acridine_pure_bulk_melt1',
    'min_lattice_length': 40,
    'prep_crystal_in_melt': False,
    'equil_time': int(2.5e5),
    'melt_temperature': 2000,
}

'''run configs to loop over - will take all combinations of the below lists - grows combinatorially!!'''

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[10, 10, 10]],
    'temperature': [100],
    'structure_identifier': ["NICOAM16"],  # , "NICOAM17"],
    'defect_rate': [0],
    'gap_rate': [0],
    'scramble_rate': [0],
    'seed': [1],
    'damping': [str(100.0)],
    'max_sphere_radius': [20],

    # static items - DO NOT SET AS LIST
    'cluster_type': 'supercell',
    'run_time': int(1e5),
    'box_type': 'p',
    'integrator': 'nosehoover',
    'print_steps': int(1e2),
    'min_inter_cluster_distance': None,
    'bulk_crystal': True,
    'machine': 'cluster',
    'run_name': 'melt_cluster_dev_test2',
    'min_lattice_length': 60,
    'prep_crystal_in_melt': True,
    'equil_time': 1e5,
    'melt_temperature': 2000,
}

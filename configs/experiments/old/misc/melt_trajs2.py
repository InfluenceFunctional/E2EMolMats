'''run configs to loop over - will take all combinations of the below lists - grows combinatorially!!'''

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[10, 10, 10]],
    'temperature': [400],
    'structure_identifier': ["NICOAM13",
                             "NICOAM16",
                             "NICOAM17"],
    'defect_rate': [0],
    'gap_rate': [0],
    'scramble_rate': [0],
    'seed': [1],
    'damping': [str(100.0)],
    'max_sphere_radius': [10, 15, 20],

    # static items - DO NOT SET AS LIST
    'cluster_type': 'supercell',
    'run_time': int(1e7),
    'box_type': 'p',
    'integrator': 'nosehoover',
    'print_steps': int(1e4),
    'min_inter_cluster_distance': 500,
    'bulk_crystal': False,
    'machine': 'cluster',
    'run_name': 'melt_trajs2',
    'min_lattice_length': 40,
}

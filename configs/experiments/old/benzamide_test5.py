'''run configs to loop over - will take all combinations of the below lists - grows combinatorially!!'''

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[3, 3, 3]],
    'temperature': [250, 300, 350],
    'structure_identifier': ["NICOAM13", "NICOAM17"],
    'defect_rate': [0, 0.1, 0.2, 0.3],
    'gap_rate': [0],
    'scramble_rate': [0],
    'seed': [1, 2],
    'damping': [str(100.0)],

    # static items - DO NOT SET AS LIST
    'run_time': int(1e7),
    'box_type': 'p',
    'integrator': 'nosehoover',
    'print_steps': int(1e3),
    'min_inter_cluster_distance': 1000,
    'bulk_crystal': False,
    'machine': 'cluster',
    'run_name': 'benzamide_test5',
}
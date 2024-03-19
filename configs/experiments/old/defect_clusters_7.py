'''run configs to loop over - will take all combinations of the below lists - grows combinatorially!!'''

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[10, 10, 10]],
    'temperature': [200, 250],
    'structure_identifier': ["NICOAM13", "NICOAM16", "NICOAM17"],
    'defect_rate': [0, 0.025, 0.05, 0.075, 0.1],
    'gap_rate': [0],
    'scramble_rate': [0],
    'seed': [1, 2],
    'damping': [str(100.0)],
    'max_sphere_radius': [15, 17.5, 20, 22.5, 25],

    # static items - DO NOT SET AS LIST
    'cluster_type': 'spherical',
    'run_time': int(1e7),
    'box_type': 'p',
    'integrator': 'nosehoover',
    'print_steps': int(1e3),
    'min_inter_cluster_distance': 1000,
    'bulk_crystal': False,
    'machine': 'cluster',
    'run_name': 'defect_clusters_7',
    'min_lattice_length': None,
}

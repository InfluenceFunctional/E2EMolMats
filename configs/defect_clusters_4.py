'''run configs to loop over - will take all combinations of the below lists - grows combinatorially!!'''

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[10, 10, 10]],
    'temperature': [250, 300],
    'structure_identifier': ["NICOAM13","NICOAM16","NICOAM17"],
    'defect_rate': [0, 0.05, 0.1, 0.15, 0.2],
    'gap_rate': [0],
    'scramble_rate': [0],
    'seed': [1, 2, 3, 4],
    'damping': [str(100.0)],
    'max_sphere_radius': [12.5, 17.5],

    # static items - DO NOT SET AS LIST
    'cluster_type': 'spherical',
    'run_time': int(1e7),
    'box_type': 'p',
    'integrator': 'nosehoover',
    'print_steps': int(1e3),
    'min_inter_cluster_distance': 1000,
    'bulk_crystal': False,
    'machine': 'cluster',
    'run_name': 'defect_clusters_4',
}
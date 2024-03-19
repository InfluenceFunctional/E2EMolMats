'''run configs to loop over - will take all combinations of the below lists - grows combinatorially!!'''

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[10, 10, 10]],
    'temperature': [100, 400],
    'structure_identifier': ["NICOAM13", "NICOAM16"],  # , "NICOAM17"],
    'defect_rate': [0],
    'gap_rate': [0],
    'scramble_rate': [0],
    'seed': [2],
    'damping': [str(100.0)],
    'max_sphere_radius': [15, 20],

    # static items - DO NOT SET AS LIST
    'cluster_type': 'supercell',
    'run_time': int(1e6),
    'box_type': 'p',
    'integrator': 'nosehoover',
    'print_steps': int(1e3),
    'min_inter_cluster_distance': None,
    'bulk_crystal': True,
    'machine': 'cluster',
    'run_name': 'crystal_in_melt_test3',
    'min_lattice_length': 60,
    'prep_crystal_in_melt': True,
    'equil_time': int(5e4),
    'melt_temperature': 1000,
}

'''run configs to loop over - will take all combinations of the below lists - grows combinatorially!!'''

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[10, 10, 10]],
    'temperature': [900],
    'structure_identifier': ["NICOAM13"],
    'defect_rate': [0],
    'gap_rate': [0],
    'scramble_rate': [0],
    'seed': [1, 2, 3, 4],
    'damping': [str(100.0)],
    'max_sphere_radius': [None],

    # static items - DO NOT SET AS LIST
    'cluster_type': 'supercell',
    'run_time': int(1e6),
    'box_type': 'p',
    'integrator': 'npt',
    'print_steps': int(1e2),
    'min_inter_cluster_distance': None,
    'bulk_crystal': True,
    'machine': 'cluster',
    'run_name': 'amorphous_bulk1',
    'min_lattice_length': 40,
}

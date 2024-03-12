'''run configs to loop over - will take all combinations of the below lists - grows combinatorially!!'''

batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[10, 10, 10]],
    'temperature': [100, 350],
    'structure_identifier': ["NICOAM07",
                             "NICOAM08",
                             "NICOAM09",
                             "NICOAM13",
                             "NICOAM14",
                             "NICOAM15",
                             "NICOAM16",
                             "NICOAM17",
                             "NICOAM18"],
    'defect_rate': [0],
    'defect_type': ['benzamide'],
    'gap_rate': [0],
    'scramble_rate': [0],
    'seed': [1],
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
    'run_name': 'new_nic_bulk_small',
    'min_lattice_length': 20,
    'prep_crystal_in_melt': False,
    'equil_time': int(2.5e5),
    'melt_temperature': 2000,
}
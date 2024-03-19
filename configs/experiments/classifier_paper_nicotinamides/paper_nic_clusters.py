'''run configs to loop over - will take all combinations of the below lists - grows combinatorially!!'''


batch_config = {
    # loop-overable (must be a list)
    'cluster_size': [[10, 10, 10]],
    'temperature': [100, 400],
    'structure_identifier': ["NICOAM13"], #, "NICOAM16", "NICOAM17"],
    'defect_rate': [0],
    'gap_rate': [0],
    'scramble_rate': [0],
    'seed': [1],
    'damping': [str(100.0)],
    'max_sphere_radius': [20, 30],

    # static items - DO NOT SET AS LIST
    'cluster_type': 'supercell',
    'run_time': int(5e6),
    'box_type': 'p',
    'integrator': 'nosehoover',
    'print_steps': int(5e2),
    'min_inter_cluster_distance': 80,
    'bulk_crystal': False,
    'machine': 'cluster',
    'run_name': 'classifier_paper_nicotinamides',
    'min_lattice_length': 60,
    'prep_crystal_in_melt': False,
    'equil_time': int(2e5),
    'melt_temperature': 2000,
}

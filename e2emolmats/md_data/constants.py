from e2emolmats.md_data.process_molecule_data import MOLECULE_DATA

FF_PATH_DICT = {'nicotinamide': 'gaff2_nicotinamid_long.lt',
                'acridine': 'gaff2_acridine_w_additive.lt'}

MOLECULE_ATOM_TYPES_MASSES = {
    'nicotinamide': "Masses\n\n"
                    "1 1.008  # ha\n"
                    "2 1.008  # h4\n"
                    "3 1.008  # hn\n"
                    "4 14.01  # n\n"
                    "5 14.01  # nb\n"
                    "6 12.01  # c\n"
                    "7 12.01  # ca\n"
                    "8 16.00  # o\n",
    'acridine': "Masses\n\n"
                "1 12.01  # ca\n"
                "2 1.008  # ha\n"
                "3 14.01  # nb\n"
                "4 15.9994  # oh\n"
                "5 1.00794  # ho\n"
}

MOLECULE_STR_TO_TYPE = {}  # mapping from atom type string to numerical index, as defined in the masses block
for key in MOLECULE_ATOM_TYPES_MASSES.keys():
    MOLECULE_STR_TO_TYPE[key] = {}
    mass_block = MOLECULE_ATOM_TYPES_MASSES[key]
    lines = mass_block.split('\n')
    for line in lines:
        splitline = line.split(' ')
        if len(line) > 0:
            try:
                int(line[0])
                MOLECULE_STR_TO_TYPE[key][splitline[-1]] = int(splitline[0])
            except ValueError:
                pass

MOLECULE_STR_TO_ATOMIC_NUM = {
    'h': 1,
    'c': 6,
    'n': 7,
    'o': 8
}  # mapping from atom type string to numerical index, as defined in the masses block

MOLECULE_NUM_ATOM_TYPES = {key: len(values) for key, values in MOLECULE_STR_TO_TYPE.items()}

MOLECULE_NUM_ATOMS = {'acridine': 23,
                      'anthracene': len(MOLECULE_DATA['anthracene']['types']),
                      '2,7-dihydroxynaphthalene': len(MOLECULE_DATA['2,7-dihydroxynaphthalene']['types']),
                      'nicotinamide': 15,
                      'benzamide': 16,
                      'isonicotinamide': 15,
                      }

# assumes straightforward 1,2,3... indexing
MOLECULE_SYM_INDICES = {key: {ind: ind for ind in range(1, MOLECULE_NUM_ATOMS[key] + 1)} for key in
                        MOLECULE_NUM_ATOMS.keys()}

MOLECULE_SHORTHAND = {'nicotinamide': 'nic1',
                      'acridine': 'AC1'}

ATOM_TYPES = {
    'nicotinamide': {
        1: 7,
        2: 7,
        3: 5,
        4: 7,
        5: 7,
        6: 7,
        7: 6,
        8: 8,
        9: 4,
        10: 3,
        11: 3,
        12: 2,
        13: 2,
        14: 1,
        15: 1
    },
    'benzamide': {
        1: 7,
        2: 7,
        3: 7,
        4: 7,
        5: 7,
        6: 7,
        7: 6,
        8: 8,
        9: 4,
        10: 3,
        11: 3,
        12: 1,
        13: 1,
        14: 1,
        15: 1,
        16: 1
    },
    'isonicotinamide': {
        1: 7,
        2: 7,
        3: 7,
        4: 5,
        5: 7,
        6: 7,
        7: 6,
        8: 8,
        9: 4,
        10: 3,
        11: 3,
        12: 2,
        13: 2,
        14: 1,
        15: 1
    },
    'acridine': {
        1: 1,
        2: 2,
        3: 2,
        4: 2,
        5: 2,
        6: 2,
        7: 2,
        8: 2,
        9: 2,
        10: 2,
        11: 2,
        12: 2,
        13: 2,
        14: 2,
        15: 3,
        16: 3,
        17: 3,
        18: 3,
        19: 3,
        20: 3,
        21: 3,
        22: 3,
        23: 3
    },
    'anthracene': {int(key): MOLECULE_STR_TO_TYPE['acridine'][val] for key, val in
                   zip(MOLECULE_DATA['anthracene']['indices'],
                       MOLECULE_DATA['anthracene']['types'])},
    '2,7-dihydroxynaphthalene': {int(key): MOLECULE_STR_TO_TYPE['acridine'][val] for key, val in
                                 zip(MOLECULE_DATA['2,7-dihydroxynaphthalene']['indices'],
                                     MOLECULE_DATA['2,7-dihydroxynaphthalene']['types'])},
}

ATOM_CHARGES = {
    'nicotinamide': {
        1: -0.326615,
        2: 0.380568,
        3: -0.585364,
        4: 0.384166,
        5: -0.38538,
        6: 0.173788,
        7: 0.6054,
        8: -0.479634,
        9: -0.779885,
        10: 0.357505,
        11: 0.357505,
        12: 0.027993,
        13: 0.034858,
        14: 0.157212,
        15: 0.077882
    },
    'benzamide': {
        1: -0.073516,
        2: -0.081906,
        3: -0.123099,
        4: -0.104813,
        5: -0.123099,
        6: -0.081906,
        7: 0.58899,
        8: -0.47625,
        9: -0.787316,
        10: 0.350274,
        11: 0.350274,
        12: 0.102542,
        13: 0.119669,
        14: 0.117946,
        15: 0.119669,
        16: 0.102542
    },
    'isonicotinamide': {
        1: 0.228047,
        2: -0.429772,
        3: 0.429792,
        4: -0.608913,
        5: 0.429792,
        6: -0.429772,
        7: 0.560881,
        8: -0.462231,
        9: -0.809865,
        10: 0.366217,
        11: 0.366217,
        12: 0.153507,
        13: 0.026297,
        14: 0.026297,
        15: 0.153507
    },
    'acridine': {
        1: -0.693636,
        2: 0.534823,
        3: -0.002672,
        4: -0.20686,
        5: -0.002672,
        6: 0.534823,
        7: -0.274751,
        8: -0.100851,
        9: -0.140261,
        10: -0.187088,
        11: -0.187088,
        12: -0.140261,
        13: -0.100851,
        14: -0.274751,
        15: 0.151995,
        16: 0.152509,
        17: 0.128637,
        18: 0.130168,
        19: 0.133735,
        20: 0.133735,
        21: 0.130168,
        22: 0.128637,
        23: 0.152509
    },
    'anthracene': {int(key): float(val) for key, val in
                   zip(MOLECULE_DATA['anthracene']['indices'],
                       MOLECULE_DATA['anthracene']['charges'])},
    '2,7-dihydroxynaphthalene': {int(key): float(val) for key, val in
                                 zip(MOLECULE_DATA['2,7-dihydroxynaphthalene']['indices'],
                                     MOLECULE_DATA['2,7-dihydroxynaphthalene']['charges'])},
}

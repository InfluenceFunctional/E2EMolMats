import os

import numpy as np


def process_thermo_data():
    f = open('screen.log', "r")
    text = f.read()
    lines = text.split('\n')
    f.close()

    skip = True
    results_dict = {'time step': [],
                    'temp': [],
                    'E_pair': [],
                    'E_mol': [],
                    'E_tot': [],
                    'Press': [],
                    'ns_per_day': []}

    if "Total wall time" not in text:  # skip analysis if the run crashed
        for key in results_dict.keys():
            results_dict[key] = np.zeros(1)
        return results_dict

    for ind, line in enumerate(lines):
        if 'ns/day' in line:
            text = line.split('ns/day')
            ns_per_day = float(text[0].split(' ')[1])
            results_dict['ns_per_day'] = ns_per_day

        if skip:
            if "Step" in line:
                skip = False
                # print(ind)
        else:
            if "Loop" in line:
                skip = True
                # print(ind)

            if not skip:
                split_line = line.split(' ')
                entries = [float(entry) for entry in split_line if entry != '']
                for ind2, key in enumerate(results_dict.keys()):
                    if key != 'ns_per_day':
                        results_dict[key].append(entries[ind2])

    for key in results_dict.keys():
        results_dict[key] = np.asarray(results_dict[key])

    if os.path.exists('tmp.out'):  # molecule-wise temperature analysis

        f = open('tmp.out', "r")
        text = f.read()
        lines = text.split('\n')
        f.close()

        frames = {}
        frame_data = []  # temp kecom internal
        for ind, line in enumerate(lines):
            if line == '\n':
                pass
            elif len(line.split()) == 0:
                pass
            elif line[0] == '#':
                pass
            elif len(line.split()) == 2:
                if len(frame_data) > 0:
                    frames[frame_num] = frame_data
                a, b = line.split()
                frame_num = int(a)
                n_mols = int(b)
                frame_data = np.zeros((n_mols, 3))
            else:
                mol_num, temp, kecom, internal = np.asarray(line.split()).astype(float)
                frame_data[int(mol_num) - 1] = temp, kecom, internal

        results_dict['thermo_trajectory'] = np.asarray(list(frames.values()))
        # averages over molecules
        results_dict['molwise_mean_temp'] = np.mean(results_dict['thermo_trajectory'][..., 0], axis=1)
        results_dict['molwise_mean_kecom'] = np.mean(results_dict['thermo_trajectory'][..., 1], axis=1)
        results_dict['molwise_mean_internal'] = np.mean(results_dict['thermo_trajectory'][..., 2], axis=1)

    return results_dict

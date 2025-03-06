import os

import numpy as np
import pandas as pd
import wandb

from e2emolmats.processing.utils import process_thermo_data, make_thermo_figs

atoms_per_molecule = {
    'nicotinamide': 15,
    'acridine': 23
}


def process_run(df, run_dir, skip_molwise_thermo, config, runs_dict, battery_full_path):
    print(f'Processing {run_dir}')
    if 'daisuke' in os.getcwd():  # todo deprecate this eventually - relevant for old runs only and will probably not work in future
        dir_split = run_dir.split('\\')
        seed = int(dir_split[0].split('seed')[-1])
        temp = int(dir_split[1].split('T')[-1])
        polymorph = dir_split[-1]
        run_config = {
            'seed': [seed],
            'temperature': [temp],
            'polymorph': [polymorph],
        }
    else:
        run_config = np.load('run_config.npy', allow_pickle=True).item()
    thermo_results_dict, analysis_code = process_thermo_data(
        run_config,
        skip_molwise_thermo,
        enforce_new_analysis=not config.latents_analysis and not config.lattice_energy_analysis
    )
    runs_dict[run_dir] = [analysis_code, run_config]
    if analysis_code != 'Thermo analysis succeeded':
        print(f'Processing {run_dir} failed ' + analysis_code)
        #continue
        return df  # if processing failed, skip this run

    thermo_figs_dict = make_thermo_figs(thermo_results_dict, run_config)
    if config.log_to_wandb:
        wandb.log(thermo_figs_dict)
    '''save results'''
    if 'num_atoms' in thermo_results_dict.keys():
        num_mols = int(thermo_results_dict['num_atoms'] / atoms_per_molecule[config.molecule])
    else:
        num_mols = thermo_results_dict['thermo_trajectory'].shape[1]
    new_row = {"run_num": run_dir,
               'num_molecules': [num_mols],
               'run_config': [run_config],
               }
    for key in run_config.keys():
        new_row.update({key: [run_config[key]]})
    for key in thermo_results_dict.keys():
        new_row.update({key: [thermo_results_dict[key]]})
    new_row.update({'time step': [thermo_results_dict['time step']]})
    df = pd.concat([df, pd.DataFrame.from_dict(new_row)])
    df.to_pickle(battery_full_path + '/results_df')
    return df

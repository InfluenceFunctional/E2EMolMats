"""central script for the analysis of MD trajectories"""

import glob
import os

import numpy as np
import pandas as pd
import wandb

from e2emolmats.analysis.utils import atoms_per_molecule
from e2emolmats.common.utils import dict2namespace
from e2emolmats.reporting.combined_analysis import combined_trajectory_analysis
from e2emolmats.processing.utils import mode_settings, make_thermo_figs, process_thermo_data, get_melt_progress, \
    relabel_defects
from e2emolmats.reporting.utils import runs_summary_table

traj_thermo_keys = ['temp', 'E_pair',
                    'E_mol', 'E_tot', 'PotEng',
                    'Press', 'Volume']

'paths for analysis of acridine melt point'
acridine_melt_paths = [
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface1/',
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface2/', # failed
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface3/', # failed
    r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface5/',
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface6/', # something really weird happened here
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan1/',

    # old acridine ff
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface14/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface15/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface16_1/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface16_2/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface16_3/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface16_4/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface17_1/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface17_3/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface17_4/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_melt_interface18/',
    #r'D:\crystal_datasets\acridine_melt_interface19/', # anthracene
    #r'D:\crystal_datasets\acridine_melt_interface20/'  # 2,7-DHN
]
acridine_scan_paths = [
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan2/',
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan3/', # first successful scan batch, with some refreezing
    #r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan4/', # single test
    r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan5/',  # shorter test to compare new thermostat

]
'paths for analysis of nicotinamide melt point'
# battery_paths = [
#     r'D:\crystal_datasets\nic_melt_interface1/',
#     r'D:\crystal_datasets\nic_melt_interface2/',
#     r'D:\crystal_datasets\nic_melt_interface3/'
# ]
'paths for analysis of nicotinamide cluster stability'
# battery_paths = [
#     r'D:\crystal_datasets\nic_cluster1/',
#     r'D:\crystal_datasets\nic_cluster2/'
# ]
'paths for analysis of acridine cluster stability'
acridine_cluster_paths = [
    # old acridine ff
    # r'D:\crystal_datasets\acridine_cluster4/',
    # r'D:\crystal_datasets\acridine_cluster5/',
    # r'D:\crystal_datasets\acridine_cluster6/',
    # r'D:\crystal_datasets\acridine_cluster7/',
    # r'D:\crystal_datasets\acridine_cluster8/',
    # r'D:\crystal_datasets\acridine_cluster9/',
    # r'D:\crystal_datasets\acridine_cluster10/',
    # r'D:\crystal_datasets\acridine_cluster11/',
    # r'D:\crystal_datasets\acridine_cluster12/',
    # r'D:\crystal_datasets\acridine_cluster13/',  # form 9 melt fix
    # r'D:\crystal_datasets\acridine_cluster14/',  # long runs
    # r'D:\crystal_datasets\acridine_cluster15/',  # init 27DHN runs
    # r'D:\crystal_datasets\acridine_cluster15_retest/',  # trying to rerun 15, where many runs failed

]
'paths for acridine latent heats of fusion'
acridine_latent_paths = [
    # old acridine ff
    # r'D:\crystal_datasets\acridine_latents_battery1/',
    # r'D:\crystal_datasets\acridine_latents_battery2/',
]
acridine_cp_paths = [
    'D:\crystal_datasets\daisuke_cp_runs'
]
acridine_cp2_paths = [
    # r'D:\crystal_datasets\acridine_w_new_ff\acridine_cp1',
    r'D:\crystal_datasets\acridine_w_new_ff\acridine_cp2',

    ##old acridine ff
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_cp1',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_cp2',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_cp3',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_latents_battery1/',
    # r'D:\crystal_datasets\acridine_w_old_ff/acridine_latents_battery2/',
]

acridine_lattice_energy_paths = [
    r'D:\crystal_datasets\acridine_w_new_ff\acridine_lattice_energy1',  # gas phases
    r'D:\crystal_datasets\acridine_w_new_ff\acridine_lattice_energy2',  # solids
    r'D:\crystal_datasets\acridine_w_new_ff\acridine_lattice_energy3',  # gas phases
    r'D:\crystal_datasets\acridine_w_new_ff\acridine_lattice_energy4',  # gas phases

]

MODE = 'acridine_scan'

if __name__ == '__main__':
    redo_analysis = True
    log_to_wandb = False
    skip_molwise_thermo = False

    (battery_paths, melt_scan_analysis, nanocluster_analysis,
     compute_melt_temps, latents_analysis,
     cp_analysis, cp2_analysis, lattice_energy_analysis) = mode_settings(
        MODE,
        acridine_cluster_paths,
        acridine_melt_paths,
        acridine_scan_paths,
        acridine_latent_paths,
        acridine_cp_paths,
        acridine_cp2_paths,
        acridine_lattice_energy_paths,
    )

    config_i = {
        'molecule': 'nicotinamide' if 'nic' in battery_paths[0] else 'acridine',
        'battery_paths': battery_paths,
        'redo_analysis': redo_analysis,
        'run_name': 'test_analysis',
        'compute_melt_temps': compute_melt_temps,
        'nanocluster_analysis': nanocluster_analysis,
        'latents_analysis': latents_analysis,
        'cp_analysis': cp_analysis,
        'cp2_analysis': cp2_analysis,
        'log_to_wandb': log_to_wandb,
        'lattice_energy_analysis': lattice_energy_analysis,
        'melt_scan_analysis': melt_scan_analysis,
    }
    config = dict2namespace(config_i)

    if config.log_to_wandb:
        wandb.init(config=config_i, project="E2EMolMats",
                   entity="mkilgour", tags=battery_paths,
                   )
        wandb.run.name = config.run_name
        wandb.run.save()

    combined_df = pd.DataFrame()
    runs_dict = {}
    for battery_path in battery_paths:
        'process and collect results battery-wise'
        print(battery_path)
        os.chdir(battery_path)
        battery_full_path = os.getcwd()

        if os.path.exists('results_df') and not config.redo_analysis:
            results_df = pd.read_pickle('results_df')
        else:
            results_df = pd.DataFrame(columns=['run_num'] + traj_thermo_keys)

        'get any directories or subdirectories down 3 levels'
        dirs = os.listdir() + glob.glob('*/*') + glob.glob('*/*/*')
        for run_dir in dirs:
            os.chdir(battery_full_path)
            try:
                os.chdir(run_dir)
            except NotADirectoryError:
                continue

            if os.path.exists('log.lammps'):
                skip_dir = False
            else:
                skip_dir = True

            if not skip_dir:
                if (run_dir not in results_df["run_num"].values) or config.redo_analysis:  #
                    print(f'Processing {run_dir}')
                    if 'daisuke' in os.getcwd():
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

                    'always do thermo analysis'
                    thermo_results_dict, analysis_code = process_thermo_data(
                        run_config,
                        skip_molwise_thermo
                    )
                    runs_dict[run_dir] = [analysis_code, run_config]
                    if analysis_code != 'Thermo analysis succeeded':
                        print(f'Processing {run_dir} failed ' + analysis_code)
                        continue  # if processing failed, skip this run

                    thermo_figs_dict = make_thermo_figs(thermo_results_dict, run_config)
                    if config.log_to_wandb:
                        wandb.log(thermo_figs_dict)

                    '''save results'''
                    new_row = {"run_num": run_dir,
                               'num_molecules': [thermo_results_dict['thermo_trajectory'].shape[
                                                     1]] if 'thermo_trajectory' in thermo_results_dict.keys() else int(
                                   thermo_results_dict['num_atoms'] / atoms_per_molecule[config.molecule]),
                               'run_config': [run_config],
                               }
                    for key in run_config.keys():
                        new_row.update({key: [run_config[key]]})
                    for key in thermo_results_dict.keys():
                        new_row.update({key: [thermo_results_dict[key]]})
                    new_row.update({'time step': [thermo_results_dict['time step']]})
                    results_df = pd.concat([results_df, pd.DataFrame.from_dict(new_row)])
                    results_df.to_pickle(battery_full_path + '/results_df')

        # visualize something with runs dict, maybe as a table
        # only for unfinished runs or when reprocessing
        if len(runs_dict) > 0 and config.redo_analysis:
            summary_fig = runs_summary_table(runs_dict, battery_path)
            summary_fig.show(renderer='browser')

        if config.compute_melt_temps:
            results_df.reset_index(drop=True, inplace=True)
            results_df['melt_slope'], results_df['melt_magnitude'] = get_melt_progress(results_df)
            results_df.to_pickle(battery_full_path + '/results_df')

        # collate dataframes
        if len(combined_df) > 0:
            combined_df = pd.concat([combined_df, results_df])
        else:
            combined_df = results_df
        combined_df.reset_index(drop=True, inplace=True)

    combined_df = relabel_defects(combined_df)
    combined_trajectory_analysis(config, combined_df, wandb)

    if config.log_to_wandb:
        wandb.finish()

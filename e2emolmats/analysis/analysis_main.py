"""central script for the analysis of MD trajectories"""

import glob
import os

import numpy as np
import pandas as pd
import wandb
from scipy.ndimage import gaussian_filter1d

from e2emolmats.common.utils import dict2namespace
from e2emolmats.reporting.utils import (process_thermo_data, make_thermo_fig,
                                        get_melt_progress, compute_and_plot_melt_slopes,
                                        plot_melt_points, POLYMORPH_MELT_POINTS, runs_summary_table,
                                        crystal_stability_analysis, latent_heat_analysis, analyze_heat_capacity,
                                        confirm_melt, df_row_melted, cp_and_latent_analysis)

traj_thermo_keys = ['temp', 'E_pair', 'E_mol', 'E_tot', 'PotEng',
                    'Press', 'Volume', 'molwise_mean_temp',
                    'molwise_mean_kecom', 'molwise_mean_internal']

'paths for analysis of acridine melt point'
acridine_melt_paths = [
    r'D:\crystal_datasets\acridine_melt_interface14/',
    r'D:\crystal_datasets\acridine_melt_interface15/',
    r'D:\crystal_datasets\acridine_melt_interface16_1/',
    r'D:\crystal_datasets\acridine_melt_interface16_2/',
    r'D:\crystal_datasets\acridine_melt_interface16_3/',
    r'D:\crystal_datasets\acridine_melt_interface16_4/',
    r'D:\crystal_datasets\acridine_melt_interface17_1/',
    r'D:\crystal_datasets\acridine_melt_interface17_3/',
    r'D:\crystal_datasets\acridine_melt_interface17_4/',
    r'D:\crystal_datasets\acridine_melt_interface18/',
    r'D:\crystal_datasets\acridine_melt_interface19/',
    r'D:\crystal_datasets\acridine_melt_interface20/'
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
    r'D:\crystal_datasets\acridine_cluster4/',
    r'D:\crystal_datasets\acridine_cluster5/',
    r'D:\crystal_datasets\acridine_cluster6/',
    r'D:\crystal_datasets\acridine_cluster7/',
    r'D:\crystal_datasets\acridine_cluster8/',
    r'D:\crystal_datasets\acridine_cluster9/',
    r'D:\crystal_datasets\acridine_cluster10/',
    r'D:\crystal_datasets\acridine_cluster11/',
    r'D:\crystal_datasets\acridine_cluster12/'
]
'paths for acridine latent heats of fusion'
acridine_latent_paths = [
    r'D:\crystal_datasets\acridine_latents_battery1/',
    r'D:\crystal_datasets\acridine_latents_battery2/',
]
acridine_cp_paths = [
    'D:\crystal_datasets\daisuke_cp_runs'
]
acridine_cp2_paths = [
    r'D:\crystal_datasets\acridine_cp1'
]

atoms_per_molecule = {
    'nicotinamide': 15,
    'acridine': 23
}

MODE = 'acridine_melt'

if __name__ == '__main__':
    redo_analysis = False
    log_to_wandb = False

    compute_melt_temps = False
    nanocluster_analysis = False
    latents_analysis = False
    cp_analysis = False
    cp2_analysis = False

    if MODE == 'acridine_cluster':
        battery_paths = acridine_cluster_paths
        nanocluster_analysis = True

    elif MODE == 'acridine_melt':
        battery_paths = acridine_melt_paths
        compute_melt_temps = True

    elif MODE == 'acridine_latent':
        battery_paths = acridine_latent_paths
        latents_analysis = True

    elif MODE == 'acridine_cp':
        battery_paths = acridine_cp_paths
        cp_analysis = True
        log_to_wandb = False

    elif MODE == 'acridine_cp2':
        battery_paths = acridine_cp2_paths
        cp2_analysis = True
        log_to_wandb = True

    else:
        assert False, "Unrecognized mode !!"

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
                    thermo_results_dict, analysis_code = process_thermo_data(run_config,
                                                                             True)  #skip_molwise_thermo='daisuke' in battery_paths[0])
                    runs_dict[run_dir] = [analysis_code, run_config]
                    if analysis_code != 'Thermo analysis succeeded':
                        print(f'Processing {run_dir} failed ' + analysis_code)
                        continue  # if processing failed, skip this run

                    thermo_telemetry_fig = make_thermo_fig(
                        thermo_results_dict, run_config)
                    if config.log_to_wandb:
                        wandb.log({'Thermo Data': thermo_telemetry_fig})

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
        if len(runs_dict) > 0:
            summary_fig = runs_summary_table(runs_dict, battery_path)
            summary_fig.show(renderer='browser')

        if config.compute_melt_temps:
            results_df.reset_index(drop=True, inplace=True)
            results_df['melt_slope'], results_df['melt_magnitude'] = get_melt_progress(results_df)
            results_df.to_pickle(battery_full_path + '/results_df')

        if len(combined_df) > 0:
            combined_df = pd.concat([combined_df, results_df])
        else:
            combined_df = results_df

        combined_df.reset_index(drop=True, inplace=True)

    'multi-battery analysis'
    if config.compute_melt_temps:
        combined_df = confirm_melt(combined_df)
        combined_df.drop(index=np.argwhere(combined_df['Melt Succeeded'] != True).flatten(), inplace=True)
        combined_df.reset_index(drop=True, inplace=True)

        fig, melt_estimate_dict = compute_and_plot_melt_slopes(combined_df)
        fig2 = plot_melt_points(melt_estimate_dict, POLYMORPH_MELT_POINTS[config.molecule])

        if config.log_to_wandb:
            wandb.log({'Melt Slopes': fig,
                       'Melt Temps': fig2,
                       })

    if config.nanocluster_analysis:
        combined_df = confirm_melt(combined_df)
        combined_df.drop(index=np.argwhere(combined_df['Melt Succeeded'] != True).flatten(), inplace=True)
        combined_df.reset_index(drop=True, inplace=True)
        fig, com_table, summary_fig = crystal_stability_analysis(combined_df)

        if config.log_to_wandb:
            wandb.log({'Crystal Nucleus Stability': fig,
                       'Crystal Stability Summary': com_table,
                       'Crystal Stability Averages': summary_fig})

    if config.latents_analysis:
        fig = latent_heat_analysis(combined_df)

        if config.log_to_wandb:
            wandb.log({'Latent Heat Estimation': fig})

    if config.cp_analysis:
        combined_df = analyze_heat_capacity(combined_df, atoms_per_molecule)
        aa = 1

    if config.cp2_analysis:
        fig = cp_and_latent_analysis(combined_df)

        if config.log_to_wandb:
            wandb.log({'Enthalpy Fitting': fig})


    if config.log_to_wandb:
        wandb.finish()

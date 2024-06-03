"""central script for the analysis of MD trajectories"""

import glob
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import wandb
from scipy.stats import linregress
from plotly.colors import n_colors


from e2emolmats.common.utils import dict2namespace
from e2emolmats.reporting.utils import (process_thermo_data, make_thermo_fig,
                                        get_melt_progress, compute_and_plot_melt_slopes,
                                        plot_melt_points, POLYMORPH_MELT_POINTS, runs_summary_table,
                                        crystal_stability_analysis)

traj_thermo_keys = ['temp', 'E_pair', 'E_mol', 'E_tot', 'PotEng',
                    'Press', 'Volume', 'molwise_mean_temp',
                    'molwise_mean_kecom', 'molwise_mean_internal']

'paths for analysis of acridine melt point'
# battery_paths = [
#     # r'D:\crystal_datasets\acridine_melt_interface14/',
#     # r'D:\crystal_datasets\acridine_melt_interface15/',
#     # r'D:\crystal_datasets\acridine_melt_interface16_1/',
#     # r'D:\crystal_datasets\acridine_melt_interface16_2/',
#     # r'D:\crystal_datasets\acridine_melt_interface16_3/',
#     # r'D:\crystal_datasets\acridine_melt_interface16_4/',
#     # r'D:\crystal_datasets\acridine_melt_interface17_1/',
#     # r'D:\crystal_datasets\acridine_melt_interface17_3/',
#     # r'D:\crystal_datasets\acridine_melt_interface17_4/',
#     r'D:\crystal_datasets\acridine_melt_interface18/'
# ]
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
battery_paths = [
    #r'D:\crystal_datasets\acridine_cluster1/',
    #r'D:\crystal_datasets\acridine_cluster2/',
    r'D:\crystal_datasets\acridine_cluster4/',
    r'D:\crystal_datasets\acridine_cluster5/',
    r'D:\crystal_datasets\acridine_cluster6/'

]
if __name__ == '__main__':
    config_i = {
        'molecule': 'nicotinamide' if 'nic' in battery_paths[0] else 'acridine',
        'battery_paths': battery_paths,
        'redo_analysis': False,
        'run_name': 'test_analysis',
        'log_figs': True,
        'compute_melt_temps': False,
        'CoM_analysis': True,
    }
    config = dict2namespace(config_i)

    wandb.init(config=config_i, project="E2EMolMats",
               entity="mkilgour", tags=battery_paths,
               )
    wandb.run.name = config.run_name
    wandb.run.save()

    combined_df = pd.DataFrame()
    runs_dict = {}
    for battery_path in battery_paths:
        'process and collect results battery-wise'

        os.chdir(battery_path)
        battery_full_path = os.getcwd()

        if os.path.exists('results_df') and not config.redo_analysis:
            results_df = pd.read_pickle('results_df')
        else:
            results_df = pd.DataFrame(columns=['run_num'] + traj_thermo_keys)

        'get any directories or subdirectories'
        dirs = os.listdir() + glob.glob('*/*')
        for run_dir in dirs:
            os.chdir(battery_full_path)
            try:
                os.chdir(run_dir)
            except NotADirectoryError:
                continue
            try:
                int(run_dir)
                skip_dir = False  # run directories are integer-listed
            except ValueError:
                skip_dir = True
            if not skip_dir:
                if (run_dir not in results_df["run_num"].values) or config.redo_analysis:  #
                    print(f'Processing {run_dir}')
                    run_config = np.load('run_config.npy', allow_pickle=True).item()

                    'always do thermo analysis'
                    thermo_results_dict, analysis_code = process_thermo_data(run_config)
                    runs_dict[run_dir] = [analysis_code, run_config]
                    if analysis_code != 'Thermo analysis succeeded':
                        print(f'Processing {run_dir} failed ' + analysis_code)
                        continue  # if processing failed, skip this run

                    if config.log_figs:
                        thermo_telemetry_fig, num_molecules = make_thermo_fig(traj_thermo_keys, thermo_results_dict,
                                                                              run_config)
                        wandb.log({'Thermo Data': thermo_telemetry_fig})

                    '''save results'''
                    new_row = {"run_num": run_dir,
                               "gap_rate": [run_config['gap_rate']],
                               'min_lattice_length': [run_config['min_lattice_length']],
                               'num_molecules': [thermo_results_dict['thermo_trajectory'].shape[1]],
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
        if len(runs_dict) > 0:
            summary_fig = runs_summary_table(runs_dict)
            summary_fig.show(renderer='browser')

        if config.compute_melt_temps:
            results_df.reset_index(drop=True, inplace=True)
            results_df['melt_slope'], results_df['melt_magnitude'] = get_melt_progress(results_df)
            results_df.to_pickle(battery_full_path + '/results_df')

        if len(combined_df) > 0:
            combined_df = pd.concat([combined_df, results_df])
        else:
            combined_df = results_df

    'multi-battery analysis'
    if config.compute_melt_temps:
        fig, melt_estimate_dict = compute_and_plot_melt_slopes(combined_df)
        fig2 = plot_melt_points(melt_estimate_dict, POLYMORPH_MELT_POINTS[config.molecule])

        wandb.log({'Melt Slopes': fig,
                   'Melt Temps': fig2,
                   })

    if config.CoM_analysis:
        fig = crystal_stability_analysis(combined_df)

        wandb.log({'Crystal Nucleus Stability': fig})

    wandb.finish()

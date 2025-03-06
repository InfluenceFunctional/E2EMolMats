"""central script for the analysis of MD trajectories"""

import glob
import os

import pandas as pd
import wandb

from e2emolmats.analysis.results_paths import acridine_melt_paths, acridine_scan_paths, acridine_cluster_paths, \
    acridine_latent_paths, acridine_cp_paths, acridine_cp2_paths, acridine_lattice_energy_paths
from e2emolmats.analysis.utils import process_run
from e2emolmats.common.utils import dict2namespace
from e2emolmats.reporting.combined_analysis import combined_trajectory_analysis
from e2emolmats.processing.utils import mode_settings, get_melt_progress, \
    relabel_defects
from e2emolmats.reporting.utils import runs_summary_table

traj_thermo_keys = ['temp', 'E_pair', 'E_mol', 'E_tot', 'PotEng', 'Press', 'Volume']
MODE = 'acridine_melt'
"""modes
acridine_cluster    
acridine_melt
acridine_scan
acridine_latent
acridine_cp
acridine_cp2
acridine_lattice_energy
"""

if __name__ == '__main__':
    redo_analysis = False
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
        'molecule': 'nicotinamide' if 'nic' in battery_paths[0] else 'acridine',  # todo clarify this, maybe a standalone config?
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

    combined_df = pd.DataFrame()  # dataframe containing information from a batch of MD runs
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
                if (run_dir not in results_df["run_num"].values) or config.redo_analysis:  #
                    results_df = process_run(results_df, run_dir, skip_molwise_thermo, config, runs_dict, battery_full_path)

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

    try:
        combined_df = relabel_defects(combined_df)
    except KeyError:
        print("defect relabelling failed - watch out!!")
    combined_trajectory_analysis(config, combined_df, wandb)

    if config.log_to_wandb:
        wandb.finish()

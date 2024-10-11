"""central script for the analysis of MD trajectories"""

import glob
import os

import numpy as np
import pandas as pd
import wandb
from scipy.spatial.distance import cdist

from e2emolmats.common.utils import dict2namespace
from e2emolmats.reporting.utils import (process_thermo_data, make_thermo_fig,
                                        get_melt_progress, compute_and_plot_melt_slopes,
                                        plot_melt_points, POLYMORPH_MELT_POINTS, crystal_stability_analysis,
                                        latent_heat_analysis, analyze_heat_capacity,
                                        confirm_melt, cp_and_latent_analysis, relabel_defects,
                                        get_run_potential, temperature_profile_fig, com_deviation_fig,
                                        compute_and_plot_melt_slopes_com, lattice_energy_figs, runs_summary_table,
                                        mean_temp_anomaly_fig)
import plotly.graph_objects as go
from scipy.ndimage import gaussian_filter1d

traj_thermo_keys = ['temp', 'E_pair', 'E_mol', 'E_tot', 'PotEng',
                    'Press', 'Volume', 'molwise_mean_temp',
                    'molwise_mean_kecom', 'molwise_mean_internal']

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
    r'D:\crystal_datasets\acridine_w_new_ff/acridine_interface_scan2/',
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

atoms_per_molecule = {
    'nicotinamide': 15,
    'acridine': 23
}

MODE = 'acridine_scan'


def mode_settings():
    global battery_paths, melt_scan_analysis, nanocluster_analysis, \
        compute_melt_temps, skip_molwise_thermo, latents_analysis, \
        cp_analysis, log_to_wandb, cp2_analysis, lattice_energy_analysis
    if MODE == 'acridine_cluster':
        battery_paths = acridine_cluster_paths
        nanocluster_analysis = True

    elif MODE == 'acridine_melt':
        battery_paths = acridine_melt_paths
        compute_melt_temps = True
        skip_molwise_thermo = False

    elif MODE == 'acridine_scan':
        battery_paths = acridine_scan_paths
        melt_scan_analysis = True
        skip_molwise_thermo = False

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

    elif MODE == 'acridine_lattice_energy':
        battery_paths = acridine_lattice_energy_paths
        lattice_energy_analysis = True
        log_to_wandb = True
    else:
        assert False, "Unrecognized mode !!"


if __name__ == '__main__':
    redo_analysis = True
    log_to_wandb = False

    skip_molwise_thermo = True

    compute_melt_temps = False
    nanocluster_analysis = False
    latents_analysis = False
    cp_analysis = False
    cp2_analysis = False
    lattice_energy_analysis = False
    melt_scan_analysis = True

    mode_settings()

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
                    thermo_results_dict, analysis_code = process_thermo_data(run_config, skip_molwise_thermo)
                    runs_dict[run_dir] = [analysis_code, run_config]
                    if analysis_code != 'Thermo analysis succeeded':
                        print(f'Processing {run_dir} failed ' + analysis_code)
                        continue  # if processing failed, skip this run

                    thermo_telemetry_fig = make_thermo_fig(thermo_results_dict, run_config)
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
        # if len(runs_dict) > 0:
        #     summary_fig = runs_summary_table(runs_dict, battery_path)
        #     summary_fig.show(renderer='browser')

        if config.compute_melt_temps:
            results_df.reset_index(drop=True, inplace=True)
            results_df['melt_slope'], results_df['melt_magnitude'] = get_melt_progress(results_df)
            results_df.to_pickle(battery_full_path + '/results_df')

        if len(combined_df) > 0:
            combined_df = pd.concat([combined_df, results_df])
        else:
            combined_df = results_df

        combined_df.reset_index(drop=True, inplace=True)

    combined_df = relabel_defects(combined_df)

    'multi-battery analysis'
    if config.compute_melt_temps:
        #
        # import plotly.graph_objects as go
        # from scipy.ndimage import gaussian_filter1d
        #
        # fig = go.Figure()
        #
        # good_inds = np.argwhere((combined_df['structure_identifier'] == 'acridine/Form2') *
        #                         (combined_df['temperature'] == 350)).flatten()
        #
        # for indyind, ind in enumerate(good_inds):
        #     fig.add_trace(
        #         go.Scattergl(x=combined_df['time step'][ind] / 1e6,
        #                      y=gaussian_filter1d(combined_df['E_pair'][ind] / combined_df['num_molecules'][ind], 5),
        #                      marker_color='blue', name='Form2',
        #                      showlegend=indyind == 0))
        #
        # good_inds = np.argwhere((combined_df['structure_identifier'] == 'acridine/Form4') *
        #                         (combined_df['temperature'] == 350)).flatten()
        #
        # for indyind, ind in enumerate(good_inds):
        #     fig.add_trace(
        #         go.Scattergl(x=combined_df['time step'][ind] / 1e6,
        #                      y=gaussian_filter1d(combined_df['E_pair'][ind] / combined_df['num_molecules'][ind], 5),
        #                      marker_color='red', name='Form4',
        #                      showlegend=indyind == 0))
        #
        # fig.show(renderer='browser')
        combined_df = confirm_melt(combined_df)  # TODO note for run 5 the equil time was hardcoded!
        print(f"{np.sum(combined_df['Melt Succeeded'] != True)} failed to melt!")
        combined_df.drop(index=np.argwhere(combined_df['Melt Succeeded'] != True).flatten(), inplace=True)
        combined_df.reset_index(drop=True, inplace=True)

        # temperature directional profile
        dev_slopes = []
        if True:  # False:
            mean_temp_anomaly_fig(combined_df)

            for r_ind in range(len(combined_df)):
                temperature_profile_fig(combined_df, r_ind, sigma_x=1, sigma_y=2, show_fig=True)
                fig, deviation, com_dev_slope = com_deviation_fig(combined_df, r_ind, show_fig=False)
                dev_slopes.append(com_dev_slope)

        combined_df['com_deviation_slope'] = dev_slopes

        fig, melt_estimate_dict, _ = compute_and_plot_melt_slopes(combined_df)
        fig, melt_estimate_dict2 = compute_and_plot_melt_slopes_com(combined_df)

        fig2 = plot_melt_points(melt_estimate_dict, POLYMORPH_MELT_POINTS[config.molecule])
        fig3 = plot_melt_points(melt_estimate_dict2, POLYMORPH_MELT_POINTS[config.molecule])

    if config.melt_scan_analysis:
        """
        here we need something like "melting direction/progress vs T"
        """
        r_ind = 0
        row = combined_df.iloc[r_ind]
        import plotly.graph_objects as go
        from scipy.ndimage import gaussian_filter1d

        sigma = 3
        fig = go.Figure()
        fig.add_scatter(x=row['time step'],
                        y=gaussian_filter1d(row['Temp'] / row['Temp'][-1], sigma=sigma)
                        )

        fig.add_scatter(x=row['time step'],
                        y=gaussian_filter1d(row['E_pair'] / row['E_pair'][-1], sigma=sigma)
                        )

        fig.update_layout(yaxis_range=[-np.inf, 1])
        fig.show(renderer='browser')
        aa = 1

    if config.nanocluster_analysis:
        combined_df = confirm_melt(combined_df)
        combined_df.drop(index=np.argwhere(combined_df['Melt Succeeded'] != True).flatten(), inplace=True)
        combined_df.reset_index(drop=True, inplace=True)
        crystal_stability_analysis(combined_df)

    if config.latents_analysis:
        fig = latent_heat_analysis(combined_df)

    if config.cp_analysis:
        combined_df = analyze_heat_capacity(combined_df, atoms_per_molecule)

    if config.cp2_analysis:
        combined_df = confirm_melt(combined_df)
        fig = cp_and_latent_analysis(combined_df)

    if config.lattice_energy_analysis:
        lattice_energy_figs(combined_df)

        aa = 1

    if config.log_to_wandb:
        wandb.finish()

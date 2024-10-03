"""central script for the analysis of MD trajectories"""

import glob
import os

import numpy as np
import pandas as pd
import wandb

from e2emolmats.common.utils import dict2namespace
from e2emolmats.reporting.utils import (process_thermo_data, make_thermo_fig,
                                        get_melt_progress, compute_and_plot_melt_slopes,
                                        plot_melt_points, POLYMORPH_MELT_POINTS, crystal_stability_analysis,
                                        latent_heat_analysis, analyze_heat_capacity,
                                        confirm_melt, cp_and_latent_analysis, relabel_defects,
                                        get_run_potential)

traj_thermo_keys = ['temp', 'E_pair', 'E_mol', 'E_tot', 'PotEng',
                    'Press', 'Volume', 'molwise_mean_temp',
                    'molwise_mean_kecom', 'molwise_mean_internal']

'paths for analysis of acridine melt point'
acridine_melt_paths = [
    # r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface1/',
    # r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface2/',
    # r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface3/',
    r'D:\crystal_datasets\acridine_w_new_ff/acridine_melt_interface5/',

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

MODE = 'acridine_melt'

if __name__ == '__main__':
    redo_analysis = False
    log_to_wandb = False

    skip_molwise_thermo = True

    compute_melt_temps = False
    nanocluster_analysis = False
    latents_analysis = False
    cp_analysis = False
    cp2_analysis = False
    lattice_energy_analysis = False

    if MODE == 'acridine_cluster':
        battery_paths = acridine_cluster_paths
        nanocluster_analysis = True

    elif MODE == 'acridine_melt':
        battery_paths = acridine_melt_paths
        compute_melt_temps = True
        skip_molwise_thermo=False

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
        'lattice_energy_analysis': lattice_energy_analysis
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
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        from scipy.ndimage import gaussian_filter1d
        from tqdm import tqdm

        r_ind = 0
        local_temp_keys = ['Mol Temp', 'KE', 'Internal T']
        fig = make_subplots(rows=1, cols=3, subplot_titles=local_temp_keys)
        row = combined_df.iloc[r_ind]
        for p_ind in range(3):

            #    pressure_direction = ['x', 'y', 'z'].index(row['pressure_direction'])
            bins = np.linspace(np.amin(row['com_trajectory'][:, :, p_ind]),
                               np.amax(row['com_trajectory'][:, :, p_ind]),
                               100)
            num_steps = len(row['com_trajectory'])
            temp_profile = np.zeros((num_steps, len(bins), 3, 3))
            for t_ind in tqdm(range(len(row['com_trajectory']))):
                com = row['com_trajectory'][t_ind, :, p_ind]
                com_inds = np.digitize(com, bins)

                for ind in range(3):
                    local_temp = row['thermo_trajectory'][t_ind, :, ind]
                    binned_temp = np.zeros((len(bins)))
                    for b_ind in range(len(binned_temp)):
                        binned_temp[b_ind] = np.mean(local_temp[com_inds[b_ind]])
                    temp_profile[t_ind, :, ind, p_ind] = binned_temp

        fig = make_subplots(rows=3, cols=3, subplot_titles=local_temp_keys)
        for ind in range(3):
            for p_ind in range(3):
                fig.add_trace(
                    go.Heatmap(x=bins[1:], z=gaussian_filter1d(temp_profile[:, :, ind, p_ind], sigma=5, axis=1)),
                    row=p_ind + 1, col=ind + 1)
        fig.update_xaxes(title='Position (A)')
        fig.update_yaxes(title='Local Temperature')
        fig.update_layout(coloraxis_showscale=False)
        fig.show(renderer='browser')

        fig, melt_estimate_dict, melt_estimate_dict2 = compute_and_plot_melt_slopes(combined_df)
        fig2 = plot_melt_points(melt_estimate_dict, POLYMORPH_MELT_POINTS[config.molecule])
        #fig3 = plot_melt_points(melt_estimate_dict2, POLYMORPH_MELT_POINTS[config.molecule])

        if config.log_to_wandb:
            wandb.log({'Melt Slopes': fig,
                       'Melt Temps': fig2,
                       # 'Melt Temps2': fig3,
                       })

    if config.nanocluster_analysis:
        combined_df = confirm_melt(combined_df)
        combined_df.drop(index=np.argwhere(combined_df['Melt Succeeded'] != True).flatten(), inplace=True)
        combined_df.reset_index(drop=True, inplace=True)
        crystal_stability_analysis(combined_df)
        #
        # if config.log_to_wandb:
        #     wandb.log({'Crystal Nucleus Stability': fig,
        #                'Crystal Stability Summary': com_table,
        #                'Crystal Stability Averages': summary_fig})

    if config.latents_analysis:
        fig = latent_heat_analysis(combined_df)

        if config.log_to_wandb:
            wandb.log({'Latent Heat Estimation': fig})

    if config.cp_analysis:
        combined_df = analyze_heat_capacity(combined_df, atoms_per_molecule)
        aa = 1

    if config.cp2_analysis:
        combined_df = confirm_melt(combined_df)
        fig = cp_and_latent_analysis(combined_df)

        if config.log_to_wandb:
            wandb.log({'Enthalpy Fitting': fig})

    if config.lattice_energy_analysis:
        # compute energies and volumes
        mean_potentials = np.zeros(len(combined_df))
        mean_volumes = np.zeros(len(combined_df))
        for ind, row in combined_df.iterrows():
            steps = len(row['E_mol'])
            num_mols = row['num_atoms'] // 23
            mean_potentials[ind] = np.mean(row['E_mol'][-steps // 2:] + row['E_pair'][-steps // 2:]) / num_mols
            mean_volumes[ind] = np.mean(row['Volume'][-steps // 2:]) / num_mols

        combined_df['mean_potential'] = mean_potentials
        combined_df['molar_volume'] = mean_volumes

        fixed_ids = []
        for ind, row in combined_df.iterrows():
            if row['cluster_type'] == 'gas':
                fixed_ids.append('acridine')
            else:
                fixed_ids.append(row['structure_identifier'])
        combined_df['structure_identifier'] = fixed_ids

        # for each polymorph and temperature, average solid and gas phase runs,
        # and compute the lattice energy as the difference
        energy_groups = (
            combined_df.groupby(
                ['temperature', 'cluster_type', 'structure_identifier']
            )['mean_potential'].mean())
        volume_groups = (
            combined_df.groupby(
                ['temperature', 'cluster_type', 'structure_identifier']
            )['molar_volume'].mean())

        temps = np.unique(combined_df['temperature'])
        num_temps = len(temps)
        polymorphs = np.unique(combined_df['structure_identifier'])
        polymorphs = polymorphs[polymorphs != 'acridine']
        polymorphs = polymorphs[polymorphs != 'acridine/Form8']
        num_polymorphs = len(polymorphs)
        N_A = 6.023*10**23
        a3tocm3 = 10**24
        MW = 179.21726
        lattices_information = np.zeros((num_temps, num_polymorphs, 2))
        for t_ind, temp in enumerate(temps):
            # get mean gas phase energy
            gas_phase_en = energy_groups[temp]['gas']['acridine']
            for p_ind, polymorph in enumerate(polymorphs):
                # record mean molar volume
                molar_volume = volume_groups[temp]['supercell'][polymorph]
                # cubic angstroms per mol to grams per cubic cm
                density = (a3tocm3 * MW) / (molar_volume * N_A)
                lattices_information[t_ind, p_ind, 1] = density

                # get mean solid energy
                solid_phase_en = energy_groups[temp]['supercell'][polymorph]

                lattice_energy = 4.18 * (solid_phase_en - gas_phase_en)  # convert from kcals to kJ
                lattices_information[t_ind, p_ind, 0] = lattice_energy

        reference_energy = np.array([-97.425,
                            -95.3,
                            -94.45,
                            -93.175,
                            -95.25,
                            -97.65])
        reference_density = np.array([1.2837,
                             1.2611,
                             1.2313,
                             1.2369,
                             1.2644,
                             1.2684])

        import plotly.graph_objects as go
        m_colors = ['black'] + ['red' for _ in range(num_temps)]
        m_sizes = [40] + [20 for _ in range(num_temps)]
        fig = go.Figure()
        for p_ind, polymorph in enumerate(polymorphs):
            volumes = [reference_density[p_ind]] + [lattices_information[t_ind, p_ind, 1] for t_ind in range(num_temps)]
            energies = [reference_energy[p_ind]] + [lattices_information[t_ind, p_ind, 0] for t_ind in range(num_temps)]
            fig.add_scatter(
                x=volumes[:2],
                y=energies[:2],
                marker_color=m_colors[:2],
                marker_size=m_sizes[:2],
                name=polymorph,
                )
        fig.update_layout(xaxis_title='Density g/mL', yaxis_title='Potential kJ/mol')
        fig.show(renderer='browser')

        temperatures = [-10] + temps.tolist()

        fig = go.Figure()
        for p_ind, polymorph in enumerate(polymorphs):
            volumes = [reference_density[p_ind]] + [lattices_information[t_ind, p_ind, 1] for t_ind in range(num_temps)]
            energies = [reference_energy[p_ind]] + [lattices_information[t_ind, p_ind, 0] for t_ind in range(num_temps)]
            fig.add_scatter(
                x=temperatures,
                y=energies,
                marker_color=m_colors,
                marker_size=m_sizes,
                name=polymorph,
                )
        fig.update_layout(xaxis_title='Temperature K', yaxis_title='Potential kJ/mol')
        fig.show(renderer='browser')

        fig = go.Figure()
        for p_ind, polymorph in enumerate(polymorphs):
            volumes = [reference_density[p_ind]] + [lattices_information[t_ind, p_ind, 1] for t_ind in range(num_temps)]
            energies = [reference_energy[p_ind]] + [lattices_information[t_ind, p_ind, 0] for t_ind in range(num_temps)]
            fig.add_scatter(
                x=temperatures,
                y=volumes,
                marker_color=m_colors,
                marker_size=m_sizes,
                name=polymorph,
                )
        fig.update_layout(xaxis_title='Temperature K', yaxis_title='Density g/mL')
        fig.show(renderer='browser')

        aa = 1

    if config.log_to_wandb:
        wandb.finish()

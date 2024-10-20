import numpy as np
from _plotly_utils.colors import n_colors
from plotly import graph_objects as go
from plotly.subplots import make_subplots
from scipy.ndimage import gaussian_filter, gaussian_filter1d
from scipy.spatial.distance import cdist
from scipy.stats import linregress, binned_statistic_2d


def temperature_profile_fig(combined_df, r_ind, sigma_x, sigma_y, show_fig=False):
    row = combined_df.iloc[r_ind]

    if isinstance(row['thermo_trajectory'], np.ndarray):
        local_temp_keys = ['Mol Temp', 'KE', 'Internal T']
        num_x_bins = 30
        #p_ind = ['x', 'y', 'z'].index(row['pressure_direction'])
        box_len = row['min_lattice_length'] * 2
        # bins = np.linspace(-box_len, box_len * 2,  # -box_len * 0.1, box_len*1.1,
        #                    num_x_bins)
        bins = np.linspace(0, box_len,  # -box_len * 0.1, box_len*1.1,
                           num_x_bins)
        #num_steps = len(row['com_trajectory'])

        temp_profile = row['local_temperature_profile']
        fig = make_subplots(rows=1, cols=2)
        # sampling_time = row['run_config']['run_time']
        sampling_start_index = -row['run_config']['print_steps']
        prof = temp_profile[sampling_start_index:, :]
        prof[prof == 0] = np.nan
        fig.add_trace(
            go.Heatmap(x=bins[:-2] - (bins[1] - bins[0]),
                       y=row['time step'][sampling_start_index:],
                       z=gaussian_filter(prof, sigma=[sigma_y, sigma_x])),
            row=1, col=1
        )
        timepoints = np.linspace(0, len(prof) - 1, 5).astype(int)
        for t_ind in timepoints:
            fig.add_trace(
                go.Scatter(x=bins[:-2] - (bins[1] - bins[0]),
                           y=gaussian_filter(prof, sigma=[sigma_y, sigma_x])[t_ind],
                           name=f"Time={int((row['time step'][sampling_start_index:] - row['time step'][sampling_start_index])[t_ind] / 100000)}ps",
                           ),
                row=1, col=2,
            )
        temp = row['temperature']
        fig.add_shape(type='line',
                      x0=0, x1=bins[-1],
                      y0=temp, y1=temp,
                      line=dict(color='black'),
                      xref='x',
                      yref='y',
                      row=1, col=2,
                      )
        fig.update_layout(coloraxis_showscale=False,
                          xaxis1_title='Position (A)',
                          yaxis1_title='Time Step',
                          xaxis2_title='Position (A)',
                          yaxis2_title='Local T(K)',
                          title=f'{row["structure_identifier"]} at {row["temperature"]}K')
        if show_fig:
            fig.show(renderer='browser')

        '''
        #alternate version of this fig
        
        row = combined_df.iloc[r_ind]
        p_ind = ['x', 'y', 'z'].index(row['pressure_direction'])
        num_mols = len(row['com_trajectory'][0, :, 0])
        sampling_start_index = -row['run_config']['print_steps']
        fig = go.Figure()
        fig.add_scatter(x=row['com_trajectory'][sampling_start_index:, :, p_ind].flatten(),
                        y=row['time step'][sampling_start_index:].repeat(num_mols),
                        marker_color=row['thermo_trajectory'][sampling_start_index:, :, 0].flatten(),
                        mode='markers',
                        opacity=0.25,
                        marker_size=10,
                        )
        fig.show(renderer='browser')
        '''


def com_deviation_fig(combined_df, r_ind, show_fig=False):
    row = combined_df.iloc[r_ind]
    num_time_steps = len(row['com_trajectory'])
    sampling_steps = row['run_config']['print_steps']
    sampling_start_index = num_time_steps - sampling_steps
    sampling_end_index = num_time_steps
    crystal_inds = np.arange(row['melt_indices'].crystal_start_ind, row['melt_indices'].crystal_end_ind)

    deviation = com_dist_profile(row['com_trajectory'], crystal_inds, sampling_start_index, sampling_end_index)

    fig = go.Figure()
    fig.add_scatter(x=row['time step'][-sampling_steps + 6:], y=deviation[6:])
    fig.update_layout(
        xaxis_title='Time',
        yaxis_title='Com Deviation',
        title=f'{row["structure_identifier"]} at {row["temperature"]}K')
    if show_fig:
        fig.show(renderer='browser')

    lr = linregress(x=row['time step'][-sampling_steps:], y=deviation)
    return fig, deviation, lr.slope


def com_dist_profile(com_trajectory, mol_inds, sampling_start_index, sampling_stop_index):
    deviation = np.zeros(sampling_stop_index - sampling_start_index)
    init_distmat = cdist(com_trajectory[sampling_start_index, mol_inds],
                         com_trajectory[sampling_start_index, mol_inds])
    tt = 0
    for t_ind in range(sampling_start_index, sampling_stop_index):
        com_distmats = cdist(com_trajectory[t_ind, mol_inds],
                             com_trajectory[t_ind, mol_inds])
        deviation[tt] = np.mean(np.abs(com_distmats - init_distmat))

        tt += 1
    return deviation


def lattice_energy_figs(combined_df):
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
    N_A = 6.023 * 10 ** 23
    a3tocm3 = 10 ** 24
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


def mean_temp_anomaly_fig(combined_df):
    sigma_x = 1
    sigma_y = 2
    sampling_start_index = -99  # -combined_df.iloc[0]['run_config']['print_steps']
    profile0 = combined_df.iloc[0]['local_temperature_profile']

    T_clip = 50

    temp_anomaly = np.zeros((len(combined_df), -sampling_start_index, profile0.shape[-1]))
    for r_ind in range(len(combined_df)):
        row = combined_df.iloc[r_ind]
        temp_profile = np.copy(row['local_temperature_profile'])
        # sampling_time = row['run_config']['run_time']
        prof = temp_profile[sampling_start_index:, :]
        prof[prof == 0] = np.nan
        temp = row['temperature']
        temp_anomaly[r_ind] = (prof[sampling_start_index:] - temp)
    fig = make_subplots(rows=3, cols=2)
    fig.add_trace(
        go.Heatmap(
            y=row['time step'][sampling_start_index:],
            z=temp_anomaly.mean(0).clip(max=T_clip), showscale=False,
        ),
        row=1, col=1,
    )
    timepoints = np.linspace(0, (-sampling_start_index) - 1, 5).astype(int)
    for t_ind in timepoints:
        fig.add_trace(
            go.Scatter(
                y=gaussian_filter(np.nan_to_num(temp_anomaly.mean(0).clip(max=T_clip)), sigma=[sigma_x, sigma_y])[
                    t_ind],
                name=f"Time={int((row['time step'][sampling_start_index:] - row['time step'][sampling_start_index])[t_ind] / 100000)}ps",
                legendgroup=f"Time={int((row['time step'][sampling_start_index:] - row['time step'][sampling_start_index])[t_ind] / 100000)}ps",
                showlegend=True,
            ),
            row=1, col=2,
        )

    temp_anomaly = np.zeros((len(combined_df), -sampling_start_index, profile0.shape[-1]))
    for r_ind in range(len(combined_df)):
        row = combined_df.iloc[r_ind]
        temp_profile = np.copy(row['local_temperature_profile2'])
        # sampling_time = row['run_config']['run_time']
        prof = temp_profile[sampling_start_index:, :]
        prof[prof == 0] = np.nan
        temp = row['temperature']
        temp_anomaly[r_ind] = (prof[sampling_start_index:] - np.nan_to_num(prof[sampling_start_index:]).mean())
    fig.add_trace(
        go.Heatmap(
            y=row['time step'][sampling_start_index:],
            z=temp_anomaly.mean(0).clip(max=T_clip), showscale=False,
        ),
        row=2, col=1,
    )

    timepoints = np.linspace(0, (-sampling_start_index) - 1, 5).astype(int)
    for t_ind in timepoints:
        fig.add_trace(
            go.Scatter(
                y=gaussian_filter(np.nan_to_num(temp_anomaly.mean(0).clip(max=T_clip)), sigma=[sigma_x, sigma_y])[
                    t_ind],
                name=f"Time={int((row['time step'][sampling_start_index:] - row['time step'][sampling_start_index])[t_ind] / 100000)}ps",
                legendgroup=f"Time={int((row['time step'][sampling_start_index:] - row['time step'][sampling_start_index])[t_ind] / 100000)}ps",
                showlegend=False,
            ),
            row=2, col=2,
        )

    temp_anomaly = np.zeros((len(combined_df), -sampling_start_index, profile0.shape[-1]))
    for r_ind in range(len(combined_df)):
        row = combined_df.iloc[r_ind]
        temp_profile = np.copy(row['local_temperature_profile3'])
        # sampling_time = row['run_config']['run_time']
        prof = temp_profile[sampling_start_index:, :]
        prof[prof == 0] = np.nan
        temp = row['temperature']
        temp_anomaly[r_ind] = (prof[sampling_start_index:] - np.nan_to_num(prof[sampling_start_index:]).mean())
    fig.add_trace(
        go.Heatmap(
            y=row['time step'][sampling_start_index:],
            z=temp_anomaly.mean(0).clip(max=T_clip), showscale=False
        ),
        row=3, col=1,
    )

    timepoints = np.linspace(0, (-sampling_start_index) - 1, 5).astype(int)
    for t_ind in timepoints:
        fig.add_trace(
            go.Scatter(
                y=gaussian_filter(np.nan_to_num(temp_anomaly.mean(0).clip(max=T_clip)), sigma=[sigma_x, sigma_y])[
                    t_ind],
                name=f"Time={int((row['time step'][sampling_start_index:] - row['time step'][sampling_start_index])[t_ind] / 100000)}ps",
                legendgroup=f"Time={int((row['time step'][sampling_start_index:] - row['time step'][sampling_start_index])[t_ind] / 100000)}ps",
                showlegend=False,
            ),
            row=3, col=2,
        )
    fig.update_layout(coloraxis_showscale=False,
                      xaxis1_title='Position (A)',
                      yaxis1_title='Time Step',
                      title=f'Mean Temperature Profiles')
    fig.show(renderer='browser')


def extract_local_temperature_profile(
        row,
        num_x_bins,
        num_y_bins,
        temp_type_ind,
        compute_anomaly=False,
):
    box_len = row['min_lattice_length'] * 2
    p_ind = ['x', 'y', 'z'].index(row['pressure_direction'])
    com_traj = row['com_trajectory'][..., p_ind]
    com_traj -= box_len * np.floor(com_traj / box_len)
    temp_traj = row['thermo_trajectory']
    num_mols = com_traj.shape[1]
    start_time = row['time step'][-len(com_traj)]
    end_time = row['time step'][-1]
    x_bins = np.linspace(0, box_len,  # -box_len * 0.1, box_len*1.1,
                         num_x_bins + 1)
    y_bins = np.linspace(start_time, end_time,  # -box_len * 0.1, box_len*1.1,
                         num_y_bins + 1)

    if not compute_anomaly:
        temp_profile = binned_statistic_2d(
            x=com_traj.flatten(),
            y=row['time step'][-len(com_traj):, None].repeat(num_mols, 1).flatten(),
            values=temp_traj[:, :, temp_type_ind].flatten(),
            statistic='mean',
            bins=[x_bins, y_bins],
        )
    else:
        temp_profile = binned_statistic_2d(
            x=com_traj.flatten(),
            y=row['time step'][-len(com_traj):, None].repeat(num_mols, 1).flatten(),
            values=(temp_traj[:, :, temp_type_ind] - temp_traj[:, :, temp_type_ind].mean(1)[:, None]).flatten(),
            statistic='mean',
            bins=[x_bins, y_bins],
        )

    return temp_profile.statistic.T, x_bins, y_bins


def extract_local_pe_profile(
        row,
        num_x_bins,
        num_y_bins,
        temp_type_ind,
        compute_anomaly=False,
):
    box_len = row['min_lattice_length'] * 2
    p_ind = ['x', 'y', 'z'].index(row['pressure_direction'])
    com_traj = row['com_trajectory'][..., p_ind]
    com_traj -= box_len * np.floor(com_traj / box_len)
    pe_traj = row['pe_trajectory']
    num_mols = com_traj.shape[1]
    start_time = row['time step'][-len(com_traj)]
    end_time = row['time step'][-1]
    x_bins = np.linspace(0, box_len,  # -box_len * 0.1, box_len*1.1,
                         num_x_bins + 1)
    y_bins = np.linspace(start_time, end_time,  # -box_len * 0.1, box_len*1.1,
                         num_y_bins + 1)

    if not compute_anomaly:
        pe_profile = binned_statistic_2d(
            x=com_traj.flatten(),
            y=row['time step'][-len(com_traj):, None].repeat(num_mols, 1).flatten(),
            values=pe_traj[:, :, temp_type_ind].flatten(),
            statistic='mean',
            bins=[x_bins, y_bins],
        )
    else:
        pe_profile = binned_statistic_2d(
            x=com_traj.flatten(),
            y=row['time step'][-len(com_traj):, None].repeat(num_mols, 1).flatten(),
            values=(pe_traj[:, :, temp_type_ind] - pe_traj[:, :, temp_type_ind].mean(1)[:, None]).flatten(),
            statistic='mean',
            bins=[x_bins, y_bins],
        )

    return pe_profile.statistic.T, x_bins, y_bins


def extract_local_profile(
        run_config,
        time_steps,
        values_traj,
        com_traj,
        num_x_bins,
        num_y_bins,
        compute_anomaly=False,
):
    box_len = run_config['min_lattice_length'] * 2
    p_ind = ['x', 'y', 'z'].index(run_config['pressure_direction'])
    com_traj = com_traj[..., p_ind]
    com_traj -= box_len * np.floor(com_traj / box_len)
    num_mols = com_traj.shape[1]
    start_time = time_steps[-len(com_traj)]
    end_time = time_steps[-1]
    x_bins = np.linspace(0, box_len,  # -box_len * 0.1, box_len*1.1,
                         num_x_bins + 1)
    y_bins = np.linspace(start_time, end_time,  # -box_len * 0.1, box_len*1.1,
                         num_y_bins + 1)

    if not compute_anomaly:
        value_profile = binned_statistic_2d(
            x=com_traj.flatten(),
            y=time_steps[-len(com_traj):, None].repeat(num_mols, 1).flatten(),
            values=values_traj.flatten(),
            statistic='mean',
            bins=[x_bins, y_bins],
        )
    else:
        value_profile = binned_statistic_2d(
            x=com_traj.flatten(),
            y=time_steps[-len(com_traj):, None].repeat(num_mols, 1).flatten(),
            values=(values_traj - values_traj.mean(1)[:, None]).flatten(),
            statistic='mean',
            bins=[x_bins, y_bins],
        )

    return value_profile.statistic.T, x_bins, y_bins

def bulk_vs_melt_temp_fig(y_pos, nbinsx, temp_profile, temp_anomaly_profile):
    interface_ind = nbinsx // 2
    fig = make_subplots(rows=3, cols=2)
    for row_ind in range(3):
        fig.add_scatter(
            x=y_pos,
            y=temp_profile[:, :interface_ind, row_ind].mean(1),
            name='Bulk Temperature',
            legendgroup='Bulk Temperature',
            marker_color='blue',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=1
        )
        fig.add_scatter(
            x=y_pos,
            y=temp_profile[:, interface_ind:, row_ind].mean(1),
            name='Melt Temperature',
            legendgroup='Melt Temperature',
            marker_color='red',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=1
        )

        fig.add_scatter(
            x=y_pos,
            y=temp_anomaly_profile[:, :interface_ind, row_ind].mean(1),
            name='Bulk Temperature Anomaly',
            legendgroup='Bulk Temperature Anomaly',
            marker_color='blue',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=2
        )
        fig.add_scatter(
            x=y_pos,
            y=temp_anomaly_profile[:, interface_ind:, row_ind].mean(1),
            name='Melt Temperature Anomaly',
            legendgroup='Melt Temperature Anomaly',
            marker_color='red',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=2
        )
    fig.update_xaxes(title='Time (ns)')
    fig.show(renderer='browser')


def interface_temp_fig(y_pos, nbinsx, temp_profile, temp_anomaly_profile):
    interface_ind = nbinsx // 2
    mid_bulk_ind = nbinsx // 4
    mid_melt_ind = nbinsx // (3 / 2)
    fig = make_subplots(rows=3, cols=2)

    for row_ind in range(3):
        fig.add_scatter(
            x=y_pos,
            y=temp_profile[:, :interface_ind, row_ind].mean(1),
            name='Bulk Temperature',
            legendgroup='Bulk Temperature',
            marker_color='blue',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=1
        )
        fig.add_scatter(
            x=y_pos,
            y=temp_profile[:, interface_ind:, row_ind].mean(1),
            name='Melt Temperature',
            legendgroup='Melt Temperature',
            marker_color='red',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=1
        )

        fig.add_scatter(
            x=y_pos,
            y=temp_anomaly_profile[:, :interface_ind, row_ind].mean(1),
            name='Bulk Temperature Anomaly',
            legendgroup='Bulk Temperature Anomaly',
            marker_color='blue',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=2
        )
        fig.add_scatter(
            x=y_pos,
            y=temp_anomaly_profile[:, interface_ind:, row_ind].mean(1),
            name='Melt Temperature Anomaly',
            legendgroup='Melt Temperature Anomaly',
            marker_color='red',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=2
        )

    fig.update_xaxes(title='Time (ns)')
    fig.show(renderer='browser')


def all_temp_profiles(y_pos, nbinsx, temp_profile, temp_anomaly_profile, inds_to_scan):
    interface_ind = nbinsx // 2
    melt_colors = n_colors('rgb(100, 0, 0)', 'rgb(250, 100, 100)', interface_ind, colortype='rgb')
    bulk_colors = n_colors('rgb(0, 0, 100)', 'rgb(100, 100, 250)', interface_ind, colortype='rgb')

    fig = make_subplots(rows=3, cols=2,
                        subplot_titles=[
                            'Temp', 'Temp Anomaly',
                            'KE CoM', 'KE CoM Anomaly',
                            'Intra Temp', 'Intra Temp Anomaly'
                        ])

    for row_ind in range(3):
        fig.add_scatter(
            x=y_pos,
            y=temp_profile[:, :interface_ind, row_ind].mean(1),
            name='Bulk',
            legendgroup='Bulk',
            marker_color='blue',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=1
        )
        fig.add_scatter(
            x=y_pos,
            y=temp_profile[:, interface_ind:, row_ind].mean(1),
            name='Melt Temperature',
            legendgroup='Melt Temperature',
            marker_color='red',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=1
        )

        fig.add_scatter(
            x=y_pos,
            y=temp_anomaly_profile[:, :interface_ind, row_ind].mean(1),
            name='Bulk Temperature Anomaly',
            legendgroup='Bulk Temperature Anomaly',
            marker_color='blue',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=2
        )
        fig.add_scatter(
            x=y_pos,
            y=temp_anomaly_profile[:, interface_ind:, row_ind].mean(1),
            name='Melt Temperature Anomaly',
            legendgroup='Melt Temperature Anomaly',
            marker_color='red',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=2
        )

    opacity = 0.5
    for row_ind in range(3):
        for b_ind in inds_to_scan:
            fig.add_scatter(
                x=y_pos,
                y=temp_profile[:, b_ind, row_ind],
                name='Bulk Temperature',
                legendgroup='Bulk Temperature',
                marker_color=bulk_colors[b_ind],
                showlegend=False,  #row_ind == 0 and b_ind == 0,
                opacity=opacity,
                row=row_ind + 1, col=1
            )
            fig.add_scatter(
                x=y_pos,
                y=temp_profile[:, interface_ind + b_ind, row_ind],
                name='Melt Temperature',
                legendgroup='Melt Temperature',
                marker_color=melt_colors[b_ind],
                showlegend=False,  #row_ind == 0 and b_ind == 0,
                opacity=opacity,
                row=row_ind + 1, col=1
            )

            fig.add_scatter(
                x=y_pos,
                y=temp_anomaly_profile[:, b_ind, row_ind],
                name='Bulk Temperature Anomaly',
                legendgroup='Bulk Temperature Anomaly',
                marker_color=bulk_colors[b_ind],
                showlegend=False,  #row_ind == 0 and b_ind == 0,
                opacity=opacity,
                row=row_ind + 1, col=2
            )
            fig.add_scatter(
                x=y_pos,
                y=temp_anomaly_profile[:, interface_ind + b_ind, row_ind],
                name='Melt Temperature Anomaly',
                legendgroup='Melt Temperature Anomaly',
                marker_color=melt_colors[b_ind],
                showlegend=False,  #row_ind == 0 and b_ind == 0,
                opacity=opacity,
                row=row_ind + 1, col=2
            )

    fig.update_xaxes(title='Time (ns)')
    fig.show(renderer='browser')


def all_pe_profiles(y_pos, nbinsx, pe_profile, pe_anomaly_profile, inds_to_scan):
    interface_ind = nbinsx // 2
    melt_colors = n_colors('rgb(100, 0, 0)', 'rgb(250, 100, 100)', interface_ind, colortype='rgb')
    bulk_colors = n_colors('rgb(0, 0, 100)', 'rgb(100, 100, 250)', interface_ind, colortype='rgb')

    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=[
                            'Temp', 'Temp Anomaly',
                            'KE CoM', 'KE CoM Anomaly',
                            'Intra Temp', 'Intra Temp Anomaly'
                        ])

    for row_ind in range(1):
        fig.add_scatter(
            x=y_pos,
            y=pe_profile[:, :interface_ind].mean(1),
            name='Bulk',
            legendgroup='Bulk',
            marker_color='blue',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=1
        )
        fig.add_scatter(
            x=y_pos,
            y=pe_profile[:, interface_ind:].mean(1),
            name='Melt Temperature',
            legendgroup='Melt Temperature',
            marker_color='red',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=1
        )

        fig.add_scatter(
            x=y_pos,
            y=pe_anomaly_profile[:, :interface_ind].mean(1),
            name='Bulk Temperature Anomaly',
            legendgroup='Bulk Temperature Anomaly',
            marker_color='blue',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=2
        )
        fig.add_scatter(
            x=y_pos,
            y=pe_anomaly_profile[:, interface_ind:].mean(1),
            name='Melt Temperature Anomaly',
            legendgroup='Melt Temperature Anomaly',
            marker_color='red',
            showlegend=row_ind == 0,
            row=row_ind + 1, col=2
        )

    opacity = 0.5
    for row_ind in range(1):
        for b_ind in inds_to_scan:
            fig.add_scatter(
                x=y_pos,
                y=pe_profile[:, b_ind],
                name='Bulk Temperature',
                legendgroup='Bulk Temperature',
                marker_color=bulk_colors[b_ind],
                showlegend=False,  #row_ind == 0 and b_ind == 0,
                opacity=opacity,
                row=row_ind + 1, col=1
            )
            fig.add_scatter(
                x=y_pos,
                y=pe_profile[:, interface_ind + b_ind],
                name='Melt Temperature',
                legendgroup='Melt Temperature',
                marker_color=melt_colors[b_ind],
                showlegend=False,  #row_ind == 0 and b_ind == 0,
                opacity=opacity,
                row=row_ind + 1, col=1
            )

            fig.add_scatter(
                x=y_pos,
                y=pe_anomaly_profile[:, b_ind],
                name='Bulk Temperature Anomaly',
                legendgroup='Bulk Temperature Anomaly',
                marker_color=bulk_colors[b_ind],
                showlegend=False,  #row_ind == 0 and b_ind == 0,
                opacity=opacity,
                row=row_ind + 1, col=2
            )
            fig.add_scatter(
                x=y_pos,
                y=pe_anomaly_profile[:, interface_ind + b_ind],
                name='Melt Temperature Anomaly',
                legendgroup='Melt Temperature Anomaly',
                marker_color=melt_colors[b_ind],
                showlegend=False,  #row_ind == 0 and b_ind == 0,
                opacity=opacity,
                row=row_ind + 1, col=2
            )

    fig.update_xaxes(title='Time (ns)')
    fig.show(renderer='browser')


def single_run_thermo_fig(row):
    good_keys = ['box_type', 'min_inter_cluster_distance', 'run_name', 'damping',
                 'defect_rate', 'defect_type', 'gap_rate', 'max_sphere_radius', 'min_lattice_length',
                 'scramble_rate', 'seed', 'structure_identifier', 'temperature']

    start_time_index = row['sampling_start_index']
    crystal_inds = np.arange(row['melt_indices'].crystal_start_ind, row['melt_indices'].crystal_end_ind)
    melt_inds = np.arange(row['melt_indices'].melt_start_ind, row['melt_indices'].melt_end_ind)
    mobil_cutoff = 4

    row['Melt Mobile Fraction'] = np.mean(row['com_mobility'][:, melt_inds] > mobil_cutoff, axis=1)
    row['Bulk Mobile Fraction'] = np.mean(row['com_mobility'][:, crystal_inds] > mobil_cutoff, axis=1)

    time = row['time step'][start_time_index:] / 1e6
    property_names = [key[5:] for key in row.keys() if 'Bulk' in key]
    sigma = 5
    fig = make_subplots(rows=1, cols=len(property_names), shared_xaxes=True, shared_yaxes=False,
                        subplot_titles=property_names)
    for ind in range(len(property_names)):
        key = 'Bulk ' + property_names[ind]
        fig.add_scatter(
            x=time,
            y=gaussian_filter1d(row[key][start_time_index:], sigma=sigma),
            name='Bulk',
            legendgroup='Bulk',
            showlegend=ind == 0,
            marker_color='blue',
            row=1, col=ind + 1
        )

        key = 'Melt ' + property_names[ind]
        fig.add_scatter(
            x=time,
            y=gaussian_filter1d(row[key][start_time_index:], sigma=sigma),
            name='Melt',
            legendgroup='Melt',
            showlegend=ind == 0,
            marker_color='red',
            row=1, col=ind + 1
        )
    fig.update_layout(title=str({key: row[key] for key in good_keys[0:5]})
                            + '<br>' + str({key: row[key] for key in good_keys[5:10]})
                            + '<br>' + str({key: row[key] for key in good_keys[10:]}),
                      title_automargin=True,
                      title_pad_t=20,
                      margin_t=180)
    fig.update_layout(yaxis4_type='log')
    fig.update_xaxes(title='Time (ns)')
    fig.show(renderer='browser')


def ramped_melt_T_extraction(row):
    # new melt point extraction
    start_time_index = row['sampling_start_index']
    crystal_inds = np.arange(row['melt_indices'].crystal_start_ind, row['melt_indices'].crystal_end_ind)
    mobil_cutoff = 4
    melt_tolerance = 0.5
    time = row['time step'][start_time_index:] / 1e6
    temp = row['Temp'][start_time_index:]
    mobility_fraction = np.mean(row['com_mobility'][start_time_index:, crystal_inds] > mobil_cutoff, axis=1)
    mobility_slope = np.diff(mobility_fraction, prepend=np.zeros(1))
    melting_flag = gaussian_filter1d(((mobility_slope > 0) * (mobility_fraction > 0.05)).astype(float),
                                     sigma=melt_tolerance, mode='nearest') >= 0.95
    if np.sum(melting_flag) > 0:
        melt_start_index = np.min(np.argwhere(melting_flag))
        melting_temp = np.mean(temp[melt_start_index:-5:melt_start_index + 5])
    else:
        melting_temp = np.nan
    return melting_temp

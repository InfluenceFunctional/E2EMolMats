import io
import os

import numpy as np
from _plotly_utils.colors import get_colorscale
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import linregress
from scipy.spatial.distance import cdist
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from scipy.ndimage import gaussian_filter1d
import plotly.express as px
import pandas as pd
import plotly.colors as pc
from plotly.colors import n_colors


def make_thermo_fig(thermo_results_dict, run_config):
    keys_to_use = []
    for key in thermo_results_dict.keys():
        if isinstance(thermo_results_dict['Temp'], np.ndarray):
            if thermo_results_dict[key].ndim == 1:
                keys_to_use += [key]

    n_keys = len(keys_to_use)
    cols = 6

    thermo_telemetry_fig = make_subplots(rows=int(np.ceil(n_keys / cols)), cols=cols, subplot_titles=keys_to_use)
    ind = 0
    for i, key in enumerate(keys_to_use):
        row = ind // cols + 1
        col = ind % cols + 1
        thermo_telemetry_fig.add_trace(
            go.Scattergl(x=thermo_results_dict['time step'] / 1e6,
                         y=thermo_results_dict[key],
                         name=key, showlegend=False),
            row=row, col=col
        )
        thermo_telemetry_fig.add_trace(
            go.Scattergl(x=thermo_results_dict['time step'] / 1e6,
                         y=gaussian_filter1d(thermo_results_dict[key], sigma=2),
                         name=key, showlegend=False),
            row=row, col=col
        )
        ind += 1

    thermo_telemetry_fig.update_xaxes(title_text="Time (ns)")
    good_keys = ['box_type', 'min_inter_cluster_distance', 'run_name', 'damping',
                 'defect_rate', 'defect_type', 'gap_rate', 'max_sphere_radius', 'min_lattice_length',
                 'scramble_rate', 'seed', 'structure_identifier', 'temperature']
    thermo_telemetry_fig.update_layout(title=str({key: run_config[key] for key in good_keys[0:5]})
                                             + '<br>' + str({key: run_config[key] for key in good_keys[5:10]})
                                             + '<br>' + str({key: run_config[key] for key in good_keys[10:]}),
                                       title_automargin=True,
                                       title_pad_t=20,
                                       margin_t=180)
    # thermo_telemetry_fig.update_layout(title=f"{run_config['structure_identifier']}, "
    #                                          f"Gap Rate {run_config['gap_rate']:.2f}, "
    #                                          f"System Size {run_config['min_lattice_length']} "
    #                                          f"Cubic Angstrom or {thermo_results_dict['thermo_trajectory'].shape[1]} Molecules"
    #                                          f"Sampling temperature {run_config['temperature']}K"
    #                                    )
    return thermo_telemetry_fig


def process_thermo_data(run_config, skip_molwise_thermo=False):
    results_dict = {}
    try:
        f = open('screen.log', "r")
        text = f.read()
    except FileNotFoundError or io.UnsupportedOperation:
        for key in results_dict.keys():
            results_dict[key] = np.zeros(1)
        return results_dict, 'screen.log missing!'
    lines = text.split('\n')
    f.close()
    hit_minimization = False
    skip = True

    if "Total wall time" not in text:  # skip analysis if the run crashed or is unfinished
        for key in results_dict.keys():
            results_dict[key] = np.zeros(1)
        return results_dict, 'Run unfinished!'
    else:
        walltime_str = text.split('Total wall time: ')[-1].split('\n')[0]
        walltime_str = walltime_str.split(':')
        total_time = int(walltime_str[-1]) + int(walltime_str[1]) * 60 + int(walltime_str[0]) * 3600
        if total_time == 0:
            return results_dict, 'Run unfinished!'
        else:
            results_dict['total_time'] = total_time
    lines = [line for line in lines if 'warning' not in line.lower()]
    for ind, line in enumerate(lines):
        if 'ns/day' in line:
            text = line.split('ns/day')
            ns_per_day = float(text[0].split(' ')[1])
            results_dict['ns_per_day'] = ns_per_day

        if not hit_minimization:
            if 'Minimization stats' in line:
                hit_minimization = True
        elif skip:
            if "Step" in line:
                split_line = line.split(' ')
                screen_keys = [entry for entry in split_line if entry != '']
                for key in screen_keys:
                    if key not in results_dict.keys():
                        results_dict[key] = []
                skip = False
        else:
            if "Loop" in line:
                splitline = line.split(' ')
                results_dict['num_atoms'] = int(splitline[-2])
                skip = True

            if not skip:
                split_line = line.split(' ')
                entries = [float(entry) for entry in split_line if entry != '']
                # order of screen keys updated in every print block so this should be robust
                for ind2, key in enumerate(screen_keys):
                    results_dict[key].append(entries[ind2])

    results_dict['time step'] = results_dict['Step']

    for key in results_dict.keys():
        results_dict[key] = np.asarray(results_dict[key])

    if os.path.exists('com.out'):  # molecule-wise temperature analysis
        num_crystal_mols = run_config['melt_indices'].crystal_end_ind - run_config['melt_indices'].crystal_start_ind + 1
        if num_crystal_mols <= 1:
            return results_dict, "Nanocrystal contains only 0-1 molecules!"

        frames = read_lammps_thermo_traj('com.out')
        results_dict['com_trajectory'] = np.asarray(list(frames.values()))
        crystal_radius_traj = results_dict['com_trajectory'][:,
                              run_config['melt_indices'].crystal_start_ind - 1:
                              run_config['melt_indices'].crystal_end_ind - 1]

        crystal_radius = np.zeros(len(crystal_radius_traj))
        crystal_min_dist = np.zeros_like(crystal_radius)
        crystal_mean_dist = np.zeros_like(crystal_radius)
        crystal_max_dist = np.zeros_like(crystal_radius)

        for ind, step in enumerate(crystal_radius_traj):
            crystal_radius[ind] = np.mean(np.linalg.norm(step - step.mean(0), axis=1))
            dmat = cdist(step, step)
            crystal_min_dist[ind] = np.amin(dmat + np.eye(dmat.shape[-1]) * 100)
            crystal_mean_dist[ind] = np.mean(dmat)
            crystal_max_dist[ind] = np.amax(dmat)

        results_dict['crystal_radius_trajectory'] = crystal_radius
        results_dict['crystal_minimum_distance'] = crystal_min_dist
        results_dict['crystal_mean_distance'] = crystal_mean_dist
        results_dict['crystal_mean_distance'] = crystal_max_dist

    if not skip_molwise_thermo:
        if os.path.exists('tmp.out'):  # molecule-wise temperature analysis
            frames = read_lammps_thermo_traj('tmp.out')
            results_dict['thermo_trajectory'] = np.asarray(list(frames.values()))
            # averages over molecules
            results_dict['molwise_mean_temp'] = np.mean(results_dict['thermo_trajectory'][..., 0], axis=1)
            results_dict['molwise_mean_kecom'] = np.mean(results_dict['thermo_trajectory'][..., 1], axis=1)
            results_dict['molwise_mean_internal'] = np.mean(results_dict['thermo_trajectory'][..., 2], axis=1)
            return results_dict, 'Thermo analysis succeeded'

        else:

            return results_dict, 'Missing thermo trajectory!'
    else:
        return results_dict, 'Thermo analysis succeeded'


def read_lammps_thermo_traj(filename):
    f = open(filename, "r")
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

    return frames


def get_melt_point(temperature, energy_traj):
    num_points = len(energy_traj)
    normed_energy_traj = energy_traj - energy_traj.min()
    normed_energy_traj /= normed_energy_traj.max()
    linreg_result = linregress(temperature[:num_points // 5], normed_energy_traj[:num_points // 5])
    normed_energy_traj2 = normed_energy_traj - (linreg_result.slope * temperature + linreg_result.intercept)

    heaviside_max = np.mean(normed_energy_traj2[-10:])
    loss = lambda b: np.sum((normed_energy_traj2 - np.heaviside(temperature - b, heaviside_max) * heaviside_max) ** 2)

    res = minimize_scalar(loss, bounds=[temperature[0], temperature[-1]], method='bounded')
    melt_temperature = res.x
    melt_temp_index = np.argmin(np.abs(temperature - melt_temperature))
    maxstep = min(len(energy_traj), melt_temp_index + 100)
    minstep = min(len(energy_traj), melt_temp_index)
    pre_step_value = np.average(normed_energy_traj2[melt_temp_index - 100: melt_temp_index])
    post_step_value = np.average(normed_energy_traj2[minstep:maxstep])
    step_size = post_step_value - pre_step_value

    melt_fit_fig = make_subplots(rows=1, cols=2)  # visualize fit
    melt_fit_fig.add_scattergl(x=temperature, y=normed_energy_traj,
                               mode='markers', name='Normalized Energy',
                               row=1, col=1)
    melt_fit_fig.add_scattergl(x=temperature, y=np.heaviside(temperature - melt_temperature,
                                                             step_size) * step_size + linreg_result.slope * temperature + linreg_result.intercept,
                               mode='markers', name='Fit',
                               row=1, col=1)
    melt_fit_fig.add_scattergl(x=temperature, y=normed_energy_traj2,
                               mode='markers', name='Normalized Delinearized Energy',
                               row=1, col=2)
    melt_fit_fig.add_scattergl(x=temperature, y=np.heaviside(temperature - melt_temperature, step_size) * step_size,
                               mode='markers', name='Fit',
                               row=1, col=2)

    melt_fit_fig.update_yaxes(title='Intermolecular Potential')
    melt_fit_fig.update_xaxes(title='temperature')
    #melt_fit_fig.show(renderer='browser')

    return step_size, melt_temperature, melt_fit_fig


def multi_ramp_fig(results_df):
    colors = n_colors('rgb(5,120,200)', 'rgb(250,50,5)', 301, colortype='rgb')
    gaps = np.linspace(0, 0.3, 301)
    gap_dict = {int(gap * 1000): ind for ind, gap in enumerate(gaps)}
    polymorphs = [thing['structure_identifier'].split('/')[-1] for thing in results_df['run_config']]
    seen_polymorph = {polymorph: False for polymorph in polymorphs}
    temp_vs_pe_fig = make_subplots(cols=2, rows=1, subplot_titles=['E_pair', 'Volume'])
    for ind in range(len(results_df)):
        gap = results_df['gap_rate'][ind]
        temp_vs_pe_fig.add_scattergl(x=results_df['temp'][ind][3:],
                                     y=results_df['E_pair'][ind][3:],
                                     line_color=colors[gap_dict[int(gap * 1000)]],
                                     marker_color=gap,
                                     mode='markers',
                                     marker_size=5,
                                     opacity=0.5,
                                     name=polymorphs[ind],
                                     legendgroup=polymorphs[ind],
                                     showlegend=True if not seen_polymorph[polymorphs[ind]] else False,
                                     #marker_colorbar=dict(len=0.1),
                                     row=1, col=1
                                     )
        temp_vs_pe_fig.add_scattergl(x=gaussian_filter1d(results_df['temp'][ind][3:], sigma=10),
                                     y=gaussian_filter1d(results_df['E_pair'][ind][3:], 10),
                                     line_color=colors[gap_dict[int(gap * 1000)]],
                                     marker_color=gap,
                                     #mode='markers',
                                     marker_size=5,
                                     opacity=1,
                                     name=polymorphs[ind] + " smoothed",
                                     legendgroup=polymorphs[ind] + " smoothed",
                                     showlegend=True if not seen_polymorph[polymorphs[ind]] else False,
                                     row=1, col=1
                                     )
        temp_vs_pe_fig.add_scattergl(x=results_df['temp'][ind][3:],
                                     y=results_df['Volume'][ind][3:],
                                     line_color=colors[gap_dict[int(gap * 1000)]],
                                     marker_color=gap,
                                     mode='markers',
                                     marker_size=5,
                                     opacity=0.5,
                                     name=polymorphs[ind],
                                     legendgroup=polymorphs[ind],
                                     showlegend=False,  #True if not seen_polymorph[polymorphs[ind]] else False,
                                     row=1, col=2
                                     )
        temp_vs_pe_fig.add_scattergl(x=gaussian_filter1d(results_df['temp'][ind][3:], 10),
                                     y=gaussian_filter1d(results_df['Volume'][ind][3:], 10),
                                     line_color=colors[gap_dict[int(gap * 1000)]],
                                     marker_color=gap,
                                     #mode='markers',
                                     marker_size=5,
                                     opacity=1,
                                     name=polymorphs[ind] + " smoothed",
                                     legendgroup=polymorphs[ind] + " smoothed",
                                     showlegend=False,  #True if not seen_polymorph[polymorphs[ind]] else False,
                                     row=1, col=2
                                     )
        seen_polymorph[polymorphs[ind]] = True
    temp_vs_pe_fig.update_xaxes(title="Temperature (K)")

    return temp_vs_pe_fig


def make_gap_vs_tm_fig(results_df):
    polymorphs = list(np.unique([conf['structure_identifier'].split('/')[-1] for conf in results_df['run_config']]))
    num_polymorphs = len(polymorphs)
    colors = n_colors('rgb(250,50,5)', 'rgb(5,120,200)', max(2, num_polymorphs), colortype='rgb')
    melt_temps, step_sizes, polymorph_names, melt_fit_figs = [], [], [], []
    gap_vs_tm_fig = make_subplots(cols=len(polymorphs), rows=1, subplot_titles=polymorphs)
    for ind, row in results_df.iterrows():
        linreg = linregress(row['temp'][3:250], row['Volume'][3:250])
        if linreg.slope < 0:
            # if the cluster collapses from the initial condition, consider it a failed run
            step_size, melt_temp, melt_fit_fig = get_melt_point(row['temp'][3:], row['E_pair'][3:])
            step_size = -1
            melt_temp = 1  # failed to melt
        else:
            step_size, melt_temp, melt_fit_fig = get_melt_point(row['temp'][3:], row['E_pair'][3:])

        polymorph_name = row['run_config']['structure_identifier'].split('/')[1]

        polymorph_names.append(polymorph_name)
        step_sizes.append(step_size)
        melt_temps.append(melt_temp)
        melt_fit_figs.append(melt_fit_fig)
        if step_size > 0.25:
            gap_vs_tm_fig.add_scattergl(y=np.asarray(melt_temp), x=np.asarray(row['gap_rate']),
                                        opacity=1,
                                        # marker_size=np.asarray(step_size) * 50,
                                        marker_color=colors[polymorphs.index(polymorph_name)],
                                        mode='markers',
                                        name=row['run_num'],
                                        showlegend=False,
                                        # legendgroup=polymorph_name,
                                        # showlegend=True if not seen_polymorphs[polymorph_name] else False,
                                        row=1, col=polymorphs.index(polymorph_name) + 1
                                        )
    results_df['polymorph_name'] = polymorph_names
    results_df['melt_temp'] = melt_temps
    results_df['step_size'] = step_sizes
    polymorph_inds = []
    polymorph_tm = []
    for ind, polymorph in enumerate(polymorphs):
        good_inds = np.argwhere(results_df['polymorph_name'] == polymorph).flatten()
        gap_rates = np.asarray(results_df.iloc[good_inds]['gap_rate'])
        melt_temps = np.asarray(results_df.iloc[good_inds]['melt_temp'])
        step_sizes = np.asarray(results_df.iloc[good_inds]['step_size'])
        seeds = np.asarray([rc['seed'] for rc in results_df.iloc[good_inds]['run_config']])

        good = np.argwhere(step_sizes >= 0.25).flatten()
        unique_gaps = np.unique(gap_rates[good])
        mean_tms = np.zeros(len(unique_gaps))
        for gap_ind, gap_rate in enumerate(unique_gaps):
            tms = melt_temps[good][gap_rates[good] == gap_rate]
            mean_tms[gap_ind] = tms.mean()

        linear_fit = linregress(unique_gaps[-3:], mean_tms[-3:])

        xmin = 0
        xmax = max(gap_rates)
        gap_vs_tm_fig.add_scattergl(x=[xmin, xmax],
                                    y=[linear_fit.intercept + xmin, linear_fit.intercept + xmax * linear_fit.slope],
                                    showlegend=False, row=1, col=ind + 1)

        gap_vs_tm_fig.add_scattergl(x=unique_gaps,
                                    y=mean_tms,
                                    marker_size=10,
                                    mode='markers',
                                    showlegend=False, row=1, col=ind + 1)

        polymorph_inds.extend(good_inds)
        polymorph_tm.extend([mean_tms[-1]] * len(
            good_inds))  # [linear_fit.intercept] * len(good_inds))

    tm_array = np.zeros(len(results_df))
    tm_array[polymorph_inds] = polymorph_tm
    results_df['polymorph_melt_temp'] = tm_array
    gap_vs_tm_fig.update_layout(legend={'itemsizing': 'constant'})
    gap_vs_tm_fig.update_xaxes(title="Gap Rate")
    gap_vs_tm_fig.update_yaxes(title="Melt Temp(K)")

    return gap_vs_tm_fig, results_df, melt_fit_figs


def get_melt_progress(results_df):
    melt_slopes = np.zeros(len(results_df))
    melt_magnitudes = np.zeros(len(results_df))
    for ind, row in results_df.iterrows():
        equil_time = row['run_config']['equil_time']
        run_time = row['run_config']['run_time']

        crystal_reference_time = equil_time
        crystal_time_index = np.argmin(np.abs(row['time step'] - crystal_reference_time))

        melt_reference_time = 2 * equil_time
        melt_time_index = np.argmin(np.abs(row['time step'] - melt_reference_time))

        sampling_start_time = 5 * equil_time
        sampling_start_index = np.argmin(np.abs(row['time step'] - sampling_start_time))

        sampling_end_time = 5 * equil_time + run_time
        sampling_end_index = np.argmin(np.abs(row['time step'] - sampling_end_time))

        inter_energy = row['E_pair']
        crystal_energy = inter_energy[crystal_time_index]
        melt_energy = inter_energy[melt_time_index]
        sampling_energy = inter_energy[sampling_start_index:sampling_end_index]

        lr = linregress(row['time step'][sampling_start_index:sampling_end_index], sampling_energy)
        melt_slopes[ind] = lr.slope
        melt_magnitudes[ind] = (sampling_energy[-10:].mean() - crystal_energy) / (melt_energy - crystal_energy)

    return melt_slopes, melt_magnitudes


def compute_and_plot_melt_slopes(df, show_fig=True):
    polymorph_names = []
    seeds = []
    for _, row in df.iterrows():
        polymorph_names.append(row['run_config']['structure_identifier'].split('/')[1])
        seeds.append(row['run_config']['seed'])
    df['polymorph_name'] = polymorph_names
    df['seed'] = seeds
    df.reset_index(drop=True, inplace=True)
    seeds = list(np.unique(df['seed']))
    polymorphs = list(
        np.unique([conf['structure_identifier'].split('/')[-1] for conf in df['run_config']]))
    defect_types = np.unique(df['defect_type'])
    num_polymorphs = len(polymorphs)
    num_defects = len(defect_types)
    colors = px.colors.qualitative.G10
    seen_polymorph = {polymorph: False for polymorph in polymorphs}
    min_temp = np.amin([df.iloc[ind]['run_config']['temperature'] for ind in range(len(df))])
    max_temp = np.amax([df.iloc[ind]['run_config']['temperature'] for ind in range(len(df))])
    temprange = np.linspace(min_temp, max_temp, 1000)
    melt_temps = {defect_type: {} for defect_type in defect_types}
    melt_temps2 = {defect_type: {} for defect_type in defect_types}
    fig = make_subplots(rows=num_defects, cols=2,
                        subplot_titles=['Normed Intermolecular Energy', 'Intermolecular Energy Slope'])
    for polymorph in polymorphs:
        for d_ind, defect_type in enumerate(defect_types):
            row = d_ind + 1
            good_inds = np.argwhere((df['polymorph_name'] == polymorph) * (df['defect_type'] == defect_type)).flatten()
            if len(good_inds) > 0:
                temperatures = np.asarray([elem['temperature'] for elem in df.iloc[good_inds]['run_config']])
                melt_magnitudes = np.asarray(df.iloc[good_inds]['melt_magnitude']).flatten()
                melt_slopes = np.asarray(df.iloc[good_inds]['melt_slope']).flatten()

                temps = np.unique(temperatures)
                mag_at_t = np.array(
                    [np.mean([melt_magnitudes[i] for i in range(len(melt_magnitudes)) if temperatures[i] == temp]) for
                     temp
                     in
                     temps])
                slope_at_t = np.array(
                    [np.mean([melt_slopes[i] for i in range(len(melt_slopes)) if temperatures[i] == temp]) for temp in
                     temps])

                mag_spline = np.interp(temprange, temps, np.maximum.accumulate(mag_at_t))
                slope_spline = np.interp(temprange, temps, np.maximum.accumulate(slope_at_t))

                melt_T = temprange[np.argmin(np.abs(slope_spline))]
                inter_temps = (temps[0:-1] + (temps[1:] - temps[0:-1]))
                diffs = np.diff(mag_at_t) / inter_temps
                melt_T2 = inter_temps[np.argmax(diffs)]
                melt_temps[defect_type][polymorph] = melt_T
                melt_temps2[defect_type][polymorph] = melt_T2
                fig.add_scattergl(x=temperatures,
                                  y=melt_magnitudes,
                                  mode='markers',
                                  name=polymorph,
                                  legendgroup=polymorph,
                                  marker_size=7,
                                  opacity=0.5,
                                  showlegend=False,
                                  marker_color=colors[polymorphs.index(polymorph)],
                                  row=row, col=1
                                  )

                fig.add_scattergl(x=temperatures,
                                  y=melt_slopes,
                                  mode='markers',
                                  name=polymorph,
                                  legendgroup=polymorph,
                                  marker_size=7,
                                  opacity=0.5,
                                  showlegend=False,
                                  marker_color=colors[polymorphs.index(polymorph)],
                                  row=row, col=2
                                  )

                fig.add_scattergl(x=temps,
                                  y=mag_at_t,
                                  mode='markers',
                                  name=polymorph,
                                  legendgroup=polymorph,
                                  marker_size=15,
                                  showlegend=True if not seen_polymorph[polymorph] else False,
                                  marker_color=colors[polymorphs.index(polymorph)],
                                  row=row, col=1
                                  )

                fig.add_scattergl(x=temps,
                                  y=slope_at_t,
                                  mode='markers',
                                  name=polymorph,
                                  legendgroup=polymorph,
                                  marker_size=15,
                                  showlegend=False,
                                  marker_color=colors[polymorphs.index(polymorph)],
                                  row=row, col=2
                                  )

    fig.update_xaxes(title='Temperature (K)')
    if show_fig:
        fig.show(renderer='browser')
    return fig, melt_temps, melt_temps2


def plot_melt_points(melt_estimate_dict, true_melts_dict, show_fig=True):
    """
    use all data to estimate the melt point for each polymorph
    """
    from e2emolmats.analysis.free_energy_calc import experimental_polymorphs_dict
    exp_melts_dict = {key.replace(' ', ''): experimental_polymorphs_dict[key]['T_melt'] for key in
                      experimental_polymorphs_dict.keys()}
    defect_types = np.unique(list(melt_estimate_dict.keys()))
    fig2 = go.Figure()
    true_melt_temps = np.asarray(list(true_melts_dict.values()))
    fig2.add_trace(go.Bar(x=list(true_melts_dict.keys()),
                          y=list(true_melts_dict.values()),
                          text=list(true_melts_dict.values()),
                          texttemplate="%{text:.2f}",
                          name='Reference'))
    fig2.add_trace(go.Bar(x=list(exp_melts_dict.keys()),
                          y=list(exp_melts_dict.values()),
                          text=list(exp_melts_dict.values()),
                          texttemplate="%{text:.2f}",
                          name='Experiments'),)

    for defect_type in defect_types:
        melt_temps = np.asarray(list(melt_estimate_dict[defect_type].values()))
        fig2.add_trace(go.Bar(x=list(melt_estimate_dict[defect_type].keys()),
                              y=list(melt_estimate_dict[defect_type].values()),
                              text=list(melt_estimate_dict[defect_type].values()),
                              texttemplate="%{text:.2f}",
                              name=defect_type + ' Estimate'),
                       )

    fig2.update_layout(yaxis_range=[min(np.amin(melt_temps), np.amin(true_melt_temps)) - 5,
                                    max(np.amax(melt_temps), np.amax(true_melt_temps)) + 5])
    if show_fig:
        fig2.show(renderer='browser')

    return fig2


POLYMORPH_MELT_POINTS = {'acridine':
    {
        'Form2': 273 + 109.8,
        'Form3': 273 + 106.8,
        'Form4': 273 + 89,
        'Form6': 273 + 99,
        'Form7': 273 + 101,
        'Form8': 273 + 109,
        'Form9': 273 + 108.8
    },
    'nicotinamide': {"NICOAM07": 444,  # form 5
                     "NICOAM08": 381,  # form 7, Tm 381, eta
                     "NICOAM09": 377.5,  # form 8, Tm 377.5, theta
                     "NICOAM13": 402,  # form 1, Tm 402K, alpha
                     "NICOAM14": 390,  # form 2, Tm 390K, beta
                     "NICOAM15": 388,  # form 3, Tm 388K, gamma
                     "NICOAM16": 387,  # form 4, Tm 387K, delta
                     "NICOAM17": 376,  # form 9, Tm 376K, iota
                     "NICOAM18": 382.5,  # form 6, Tm 382.5K, zeta}
                     },
}


def runs_summary_table(runs_dict, run_name):
    try:
        runs = np.asarray(list(runs_dict.keys())).astype(int)
        runs = np.sort(runs)
    except ValueError:
        runs = list(runs_dict.keys())
    codes = [runs_dict[str(key)][0] for key in runs]

    if 'defect_type' in runs_dict[list(runs_dict.keys())[0]][1].keys():
        defects = [runs_dict[str(key)][1]['defect_type'] for key in runs]
        defect_rates = [runs_dict[str(key)][1]['defect_rate'] for key in runs]
        sizes = [runs_dict[str(key)][1]['max_sphere_radius'] for key in runs]
        polymorph = [runs_dict[str(key)][1]['structure_identifier'].split('/')[-1] for key in runs]
        vals = [runs, codes, defects, defect_rates, sizes, polymorph]
        keys = ['run_num', 'run_status', 'defect_type', 'defect_rate', 'cluster_radius', 'polymorph']
    else:
        polymorph = [runs_dict[str(key)][1]['polymorph'][-1] for key in runs]
        vals = [runs, codes, polymorph]
        keys = ['run_num', 'run_status', 'polymorph']

    fig = go.Figure(
        data=go.Table(header=dict(
            values=keys),
            cells=dict(values=vals)))
    fig.update_layout(title=run_name)
    return fig


def crystal_stability_analysis(combined_df):
    crystal_size = []
    stability = []
    temperature = []
    polymorphs = []
    defects = []
    defect_types = []
    for ind, row in combined_df.iterrows():
        if 'latent' in row['run_num']:
            print("skipping misplaced 'latent' run")
            continue  # skip misplaced runs

        equil_time = row['run_config']['equil_time']
        run_time = row['run_config']['run_time']
        defects.append(row['run_config']['defect_rate'])
        defect_types.append(row['defect_type'])

        polymorphs.append(row['structure_identifier'].split('/')[-1])

        sampling_start_time = 5 * equil_time
        sampling_start_index = np.argmin(np.abs(row['time step'] - sampling_start_time))

        sampling_end_time = 5 * equil_time + run_time
        sampling_end_index = np.argmin(np.abs(row['time step'] - sampling_end_time))

        crystal_size.append(row['max_sphere_radius'])
        temperature.append(row['temperature'])
        calculated_stability = (
                np.mean(row['crystal_radius_trajectory'][sampling_start_index:sampling_start_index + 5]) /
                np.amax(row['crystal_radius_trajectory'][sampling_start_index:sampling_end_index]))

        stability.append(calculated_stability)

    present_polymorphs = np.unique(polymorphs)
    present_defect_rates = np.unique(defects)
    n_columns = len(present_polymorphs)
    n_rows = len(present_defect_rates)
    predictions_list = []
    for defect_type in np.unique(defect_types):
        predictions_list.extend(
            make_critical_nucleus_size_fig(crystal_size, defect_type, defect_types, defects, n_columns,
                                           n_rows, polymorphs, present_defect_rates,
                                           present_polymorphs, stability, temperature))

    df = pd.DataFrame(predictions_list)

    scale = pc.make_colorscale(['rgb(250,200,200)', 'rgb(0,255, 150)'])
    colors = pc.sample_colorscale(scale, df['Critical Nucleus Size'] / np.amax(df['Critical Nucleus Size']))
    df['Colors'] = colors
    sum_table = go.Figure(data=[go.Table(
        header=dict(values=['Polymorph', 'Temperature', 'Defect Type', 'Defect Rate', 'Critical Nucleus Size'],
                    fill_color='paleturquoise',
                    align='left'),
        cells=dict(values=[df['Polymorph'], df['Temperature'], df['Defect Type'], df['Defect Rate'],
                           df['Critical Nucleus Size']],
                   line_color=[df['Colors']], fill_color=[df['Colors']],
                   align='left'))
    ])
    sum_table.show(renderer='browser')

    present_polymorphs = np.unique(df['Polymorph'])
    present_defect_rates = np.unique(df['Defect Rate'])
    present_defect_types = np.unique(df['Defect Type'])
    present_temperatures = np.array([320, 330, 340])
    n_polymorphs = len(present_polymorphs)
    n_defects = len(present_defect_rates)
    n_temperatures = len(present_temperatures)
    n_defect_types = len(present_defect_types)
    means_group = df.groupby(['Temperature', 'Polymorph', 'Defect Rate', 'Defect Type'])['Critical Nucleus Size'].mean()
    stability_array = np.zeros((n_polymorphs, n_defects, n_defect_types, n_temperatures))

    for p_ind, polymorph in enumerate(present_polymorphs):
        for defect_rate_index, defect_rate in enumerate(present_defect_rates):
            for defect_type_index, defect_type in enumerate(present_defect_types):
                for t_ind, temp in enumerate(present_temperatures):
                    for row_ind, row in means_group.items():
                        if int(row_ind[0]) == temp:
                            if row_ind[1] == polymorph:
                                if row_ind[2] == defect_rate:
                                    if row_ind[3] == defect_type:
                                        stability_array[p_ind, defect_rate_index, defect_type_index, t_ind] = row

    calculated_melt_temps = {
        'Form2': 395.23,
        'Form3': 394.68,  # K
        'Form4': 351.35,  # K
        'Form6': 358.82,  # K
        'Form7': 372.52,  # K
        'Form8': 371.71,  # K
        'Form9': 376.03,
    }

    stability_array[stability_array == 0] = -np.inf

    'fit overall undercooling rate for pure samples'
    crange, pure_index, undercool_func = get_undercooling_function(calculated_melt_temps,
                                                                   present_defect_types,
                                                                   present_polymorphs,
                                                                   present_temperatures,
                                                                   stability_array)

    relative_undercooled_sizes = compute_undercool_vs_defects(calculated_melt_temps,
                                                              n_defect_types,
                                                              n_defects,
                                                              n_polymorphs,
                                                              present_defect_rates,
                                                              present_defect_types,
                                                              present_polymorphs,
                                                              present_temperatures,
                                                              pure_index,
                                                              stability_array,
                                                              undercool_func)

    relative_fig = make_subplots(rows=n_defect_types - 1, cols=1)
    row = 1
    for defect_type_index in range(n_defect_types):
        if defect_type_index != pure_index:
            for defect_rate_index in range(n_defects):
                relative_fig.add_bar(x=present_polymorphs,
                                     y=relative_undercooled_sizes[defect_type_index, defect_rate_index, :],
                                     name=f'Defect {present_defect_types[defect_type_index]} Frac={present_defect_rates[defect_rate_index]:.3f}',
                                     row=row, col=1)
            row += 1
    relative_fig.update_yaxes(title='Mean Deviation from Trend')
    relative_fig.update_xaxes(title='Polymorph')
    relative_fig.show(renderer='browser')

    colors = get_colorscale(
        'mygbm')  #n_colors('rgb(0,0,255)', 'rgb(255,0,0)', len(present_polymorphs), colortype='rgb')


    for defect_type_index, defect_type in enumerate(present_defect_types):
        overall_fig = make_subplots(rows=len(present_defect_rates), cols=2,
                                    subplot_titles=['Critical Radius vs. T', 'Critical Radius vs. Undercooling'],
                                    vertical_spacing=0.07)

        for defect_rate_index, defect_rate in enumerate(present_defect_rates):
            row = defect_rate_index + 1
            if defect_rate == 0:
                defect_type_ind = pure_index
            else:
                defect_type_ind = defect_type_index * 1

            overall_fig.add_scattergl(  # marker_color=stability_array[t_ind, , :],
                y=undercool_func(crange),
                x=crange,
                name='Average Pure Trend',
                marker_color='black',
                showlegend=True if row == 1 else False,
                legendgroup='Average Pure Trend',
                row=row, col=2
            )

            for polymorph_ind, polymorph in enumerate(present_polymorphs):
                overall_fig.add_scattergl(
                    y=stability_array[polymorph_ind, defect_rate_index, defect_type_ind, :],
                    x=present_temperatures,
                    name=polymorph,
                    showlegend=True if row == 1 else False,
                    legendgroup=polymorph,
                    marker_color=colors[polymorph_ind][1],
                    row=row, col=1
                )
                undercooling_ratio = (calculated_melt_temps[polymorph] - present_temperatures) / calculated_melt_temps[
                    polymorph]
                overall_fig.add_scattergl(
                    y=stability_array[polymorph_ind, defect_rate_index, defect_type_ind, :],
                    x=undercooling_ratio,
                    name=polymorph,
                    marker_color=colors[polymorph_ind][1],
                    showlegend=False,
                    legendgroup=polymorph,
                    row=row, col=2
                )

        overall_fig.update_layout(xaxis1_title='Temperature (K)')
        overall_fig.update_layout(xaxis2_title='Undercooling Ratio')
        overall_fig.update_yaxes(title='Critical Nucleus Radius')
        overall_fig.update_layout(title=f"{defect_type} at fractions {present_defect_rates}")
        overall_fig.show(renderer='browser')

    return None


def compute_undercool_vs_defects(calculated_melt_temps, n_defect_types, n_defects, n_polymorphs,
                                 present_defect_rates, present_defect_types, present_polymorphs, present_temperatures,
                                 pure_index, stability_array, undercool_func):
    relative_undercooled_sizes = np.zeros((n_defect_types, n_defects, n_polymorphs))
    for defect_type_index, defect_type in enumerate(present_defect_types):
        for defect_rate_index, defect_rate in enumerate(present_defect_rates):
            if defect_rate == 0:
                defect_type_ind = pure_index * 1
            else:
                defect_type_ind = defect_type_index * 1

            for polymorph_ind, polymorph in enumerate(present_polymorphs):
                undercooling_ratio = (calculated_melt_temps[polymorph] - present_temperatures) / \
                                     calculated_melt_temps[polymorph]

                vals = stability_array[polymorph_ind, defect_rate_index, defect_type_ind, :]

                finite_vals = vals[np.isfinite(vals)]

                if len(finite_vals) > 0:
                    trend_sizes = undercool_func(undercooling_ratio)[np.isfinite(vals)]
                    relative_undercooled_sizes[defect_type_index, defect_rate_index, polymorph_ind] = np.mean(finite_vals - trend_sizes)

    relative_undercooled_sizes[relative_undercooled_sizes == 0] = -np.inf

    return relative_undercooled_sizes


def get_undercooling_function(calculated_melt_temps, present_defect_types, present_polymorphs, present_temperatures,
                              stability_array):
    pure_index = int(np.argwhere(present_defect_types == 'Pure').flatten())
    pure_undercools = []
    pure_sizes = []
    for polymorph_ind, polymorph in enumerate(present_polymorphs):
        stabilities = stability_array[polymorph_ind, 0, pure_index, :]
        good_inds = np.isfinite(stabilities)
        undercools = (calculated_melt_temps[polymorph] - present_temperatures) / calculated_melt_temps[polymorph]
        pure_undercools.extend(undercools[good_inds.flatten()].flatten())
        pure_sizes.extend(stabilities[good_inds.flatten()].flatten())
    sizes = np.array(pure_sizes)
    undercools = np.array(pure_undercools)
    fit = np.polyfit(undercools, sizes, 2)
    undercool_func = np.poly1d(fit)
    crange = np.linspace(min(undercools), max(undercools), 100)
    return crange, pure_index, undercool_func


def make_critical_nucleus_size_fig(crystal_size, defect_type, defect_types, defects, n_columns, n_rows, polymorphs,
                                   present_defect_rates, present_polymorphs, stability, temperature):
    fig = make_subplots(rows=n_rows, cols=n_columns, subplot_titles=present_polymorphs, horizontal_spacing=0.035,
                        vertical_spacing=0.08)
    prediction_df = []
    xspace = np.linspace(0, 100, 1000)
    colors = n_colors('rgb(0,0,255)', 'rgb(255,0,0)', len(np.unique(temperature)) + 1, colortype='rgb')
    seen_t = {temp: False for temp in np.unique(temperature)}
    for p_ind, polymorph in enumerate(present_polymorphs):
        col = p_ind + 1
        for d_ind, defect_rate in enumerate(present_defect_rates):
            row = d_ind + 1
            for t_ind, temp in enumerate(np.unique(temperature)):
                good_inds = np.where(
                    (temperature == temp) * (np.array(polymorphs) == polymorph) * (
                            np.array(defects) == defect_rate) * (np.array(defect_types) == defect_type))[0]
                stab = np.array(stability)[good_inds]
                finite_inds = np.isfinite(stab)
                stab = stab[finite_inds]
                size = np.array(crystal_size)[good_inds]
                size = size[finite_inds]

                if len(good_inds) > 1 and len(np.unique(size)) > 1:
                    # fitting function
                    fitting_func = lambda x, slope, intercept: np.minimum(1, slope * x + intercept)
                    loss = lambda x: np.sum((stab - fitting_func(size, x[0], x[1])) ** 2)

                    regress = linregress(size, stab)
                    rslope = regress.slope
                    rint = regress.intercept
                    assert np.isfinite(rslope)
                    res = minimize(loss, x0=np.array([rslope,
                                                      rint]), options={'disp': False},
                                   bounds=((None, None), (None, None)))
                    m, b = res.x

                    # plot raw scatter
                    fig.add_scattergl(x=size, y=stab, mode='markers',
                                      legendgroup=temp, name=temp, showlegend=True if not seen_t[temp] else False,
                                      marker_color=colors[t_ind],
                                      marker_size=10,
                                      opacity=0.75,
                                      row=row, col=col)
                    # plot fitting function
                    fig.add_scattergl(x=xspace, y=fitting_func(xspace, m, b), mode='lines',
                                      legendgroup=temp, name=temp, showlegend=False, marker_color=colors[t_ind],
                                      row=row, col=col)
                    if not seen_t[temp]:
                        seen_t[temp] = True

                    # extract critical size prediction
                    critical_size_prediction = xspace[np.argmin(0.975 - fitting_func(xspace, m, b))]
                    # if any clusters were observed stable below the predicted size, trust this instead
                    stable_inds = np.argwhere(stab > 0.975)
                    if len(stable_inds) > 0:
                        minimum_stable_size = np.amin(size[stable_inds])
                        if minimum_stable_size < critical_size_prediction:
                            critical_size_prediction = minimum_stable_size

                    # record and plot critical size prediction
                    prediction_df.append({'Defect Type': defect_type,
                                          'Defect Rate': defect_rate,
                                          'Polymorph': polymorph,
                                          'Temperature': temp,
                                          'Critical Nucleus Size': np.round(critical_size_prediction, 1)
                                          })
                    fig.add_scattergl(x=[critical_size_prediction], y=[1], mode='markers',
                                      legendgroup=temp, name=temp, showlegend=False, marker_color=colors[t_ind],
                                      marker_size=15, marker_line_color='black', marker_line_width=6, opacity=0.7,
                                      row=row, col=col)
                    fig.add_scattergl(x=xspace, y=[1 for _ in range(len(xspace))], showlegend=False, mode='lines',
                                      marker_color='black',
                                      row=row, col=col)
    print(prediction_df)
    fig.update_xaxes(title='Nucleus Size', range=[0, 50])
    fig.update_yaxes(title='Nucleus Stability', range=[0, 1.1], )
    fig.show(renderer='browser')
    return prediction_df


def latent_heat_analysis(combined_df):
    latents_dict = {}
    unique_idents = np.unique(combined_df['structure_identifier'])
    n_a = 6.022 * 10 ** 23
    for id, ident in enumerate(unique_idents):
        #polymorph = ident.split('/')
        good_inds = np.argwhere(combined_df['structure_identifier'] == ident).flatten()
        good_df = combined_df.iloc[good_inds]
        melted_inds = np.argwhere(good_df['prep_bulk_melt']).flatten()
        crystal_inds = np.argwhere(good_df['prep_bulk_melt'] != True).flatten()
        good_df.reset_index(drop=True, inplace=True)
        # TODO add melt confirmation
        failed_melts = []
        for run_ind, row in good_df.iterrows():
            if run_ind in melted_inds:
                if not df_row_melted(row):
                    failed_melts.append(run_ind)

        if len(failed_melts) > 0:
            good_df.drop(index=failed_melts, inplace=True)
            good_df.reset_index(drop=True, inplace=True)
            melted_inds = np.argwhere(good_df['prep_bulk_melt']).flatten()
            crystal_inds = np.argwhere(good_df['prep_bulk_melt'] != True).flatten()

        if len(melted_inds) > 0 and len(crystal_inds) > 0:
            mean_enthalpy = np.zeros(len(good_df))
            for run_ind, row in good_df.iterrows():
                if 'E_tot' in row.keys():
                    if np.sum(np.isnan(row['E_tot'])) == 0:
                        en_key = 'E_tot'
                    else:
                        en_key = 'TotEng'
                else:
                    en_key = 'TotEng'

                energy = row[en_key] * 4.184 / row['num_molecules']  # / \
                #good_df.iloc[0]['molecule_num_atoms_dict']['acridine']  # kcal/mol -> kJ/mol
                pressure = 1  #good_df.iloc[run_ind]['Press']  # atmospheres
                volume = row['Volume']
                # atm*cubic angstrom = 1.01325e-28 kJ, normed per-mol
                PV_energy = pressure * volume * n_a / row['num_molecules'] * (1.01325 * 10 ** -25) / 1000
                run_enthalpy = energy + PV_energy
                mean_enthalpy[run_ind] = np.mean(run_enthalpy[len(run_enthalpy) // 2:])

            melted_enthalpy = np.mean(mean_enthalpy[melted_inds])
            crystal_enthalpy = np.mean(mean_enthalpy[crystal_inds])

            latents_dict[ident] = melted_enthalpy - crystal_enthalpy  # in kJ per mol
        #
        # traj_thermo_keys = ['Press', 'Volume', en_key, 'PotEng', 'Temp']
        # cols = len(traj_thermo_keys)
        # thermo_telemetry_fig = make_subplots(rows=2, cols=cols, subplot_titles=traj_thermo_keys, vertical_spacing=0.075)
        # ind = 0
        # for r_ind, df_row in good_df.iterrows():
        #     for i, key in enumerate(traj_thermo_keys):
        #         col = ind % cols + 1
        #         thermo_telemetry_fig.add_trace(
        #             go.Scattergl(x=df_row['time step'] / 1e6,
        #                          y=df_row[key],
        #                          name=key, showlegend=False),
        #             row=1, col=col
        #         )
        #         thermo_telemetry_fig.add_trace(
        #             go.Scattergl(x=df_row['time step'] / 1e6,
        #                          y=gaussian_filter1d(df_row[key], sigma=10),
        #                          name=key, showlegend=False),
        #             row=2, col=col
        #         )
        #         ind += 1
        #
        # thermo_telemetry_fig.update_xaxes(title_text="Time (ns)")
        # thermo_telemetry_fig.show(renderer='browser')

    fig = go.Figure()
    fig.add_trace(
        go.Bar(x=list(latents_dict.keys()), y=list(latents_dict.values()), name='Latent Heat'))
    fig.add_trace(go.Bar(x=['Reference 383K', 'Reference2 383K'], y=[20.682, 18.58], name='Reference'))
    fig.update_layout(yaxis_title='Latent Heat of Fusion (kJ/mol)')
    fig.show(renderer='browser')

    return fig


def analyze_heat_capacity(combined_df, atoms_per_molecule):
    kb = 0.008314462618  # boltzman factor. checked(https://en.wikipedia.org/wiki/Boltzmann_constant). kJ/mol/K
    kJ = 4.184  # kJ/kcal
    # compute heat capacities
    cps = np.zeros(len(combined_df))
    cvs = np.zeros(len(combined_df))
    for ind, row in combined_df.iterrows():
        traj_len = len(row['TotEng'])
        energy_traj = row['TotEng'][traj_len // 2:] * kJ
        enthalpy_traj = row['Enthalpy'][traj_len // 2:] * kJ
        temp_traj = row['Temp'][traj_len // 2:]

        num_molecules = row['num_atoms'] / atoms_per_molecule['acridine']
        cp = (np.mean(energy_traj ** 2) - np.mean(energy_traj) ** 2) / (
                kb * np.mean(temp_traj) ** 2) / num_molecules
        cv = (np.mean(enthalpy_traj ** 2) - np.mean(enthalpy_traj) ** 2) / (
                kb * np.mean(temp_traj) ** 2) / num_molecules

        cps[ind] = cp
        cvs[ind] = cv
    combined_df['Cp'] = cps
    combined_df['Cv'] = cvs
    return combined_df


def confirm_melt(combined_df: pd.DataFrame) -> pd.DataFrame:
    """
    Check whether any given trajectory melted
    """
    'the overall volume should increase when it melts - we can check for this'
    melt_succeeded = np.zeros(len(combined_df), dtype=np.bool_)
    for row_ind, row in combined_df.iterrows():
        success = df_row_melted(row)
        melt_succeeded[row_ind] = success

    combined_df['Melt Succeeded'] = melt_succeeded
    return combined_df


def df_row_melted(row):
    vol_traj = row['Volume']
    equil_time = row['run_config']['equil_time']
    equil_time_index = np.argmin(np.abs(row['time step'] - equil_time))
    pre_heat_mean_volume = np.mean(vol_traj[1:equil_time_index])
    max_vol = np.amax(gaussian_filter1d(vol_traj, sigma=2))
    volume_ratio = np.abs(max_vol - pre_heat_mean_volume) / pre_heat_mean_volume
    success = volume_ratio > 0.05
    return success


def cp_and_latent_analysis(combined_df):
    """plot and numerically fit run enthalpy vs T"""

    enthalpies_dict = {}
    unique_idents = np.unique(combined_df['structure_identifier'])
    n_a = 6.022 * 10 ** 23
    for id, ident in enumerate(unique_idents):
        # polymorph = ident.split('/')
        good_inds = np.argwhere(
            (combined_df['structure_identifier'] == ident) * (combined_df['prep_bulk_melt'] != True)).flatten()
        good_df = combined_df.iloc[good_inds]
        good_df.reset_index(drop=True, inplace=True)

        if len(good_df) > 0:
            temps, mean_enthalpy = extract_df_enthalpies(good_df, n_a)

            enthalpies_dict[ident] = {'Temperature': temps,
                                      'Enthalpies': mean_enthalpy,
                                      }
    'for the melt fraction'
    good_inds = np.argwhere(combined_df['Melt Succeeded']).flatten()
    good_df = combined_df.iloc[good_inds]
    good_df.reset_index(drop=True, inplace=True)

    temps, mean_enthalpy = extract_df_enthalpies(good_df, n_a)

    melt_enthalpies = {'Temperature': temps,
                       'Enthalpies': mean_enthalpy,
                       }

    melts_dict = {
        'acridine/Form2': 394.68468468468467,
        'acridine/Form3': 394.68468468468467,
        'acridine/Form4': 351.3613613613614,
        'acridine/Form6': 358.8088088088088,
        'acridine/Form7': 372.5025025025025,
        'acridine/Form8': 371.7017017017017,
        'acridine/Form9': 376.02602602602605
    }

    import plotly.graph_objects as go
    from plotly.colors import n_colors

    colors = n_colors('rgb(25, 200, 255)', 'rgb(255, 128, 0)', 10, colortype='rgb')
    latents_dict = {}
    cp_coefficients_dict = {}
    cp_at_T = {}

    m_temps, m_ens = melt_enthalpies['Temperature'], melt_enthalpies['Enthalpies']
    # 'fit numerical gradient function'
    # u_temps = np.unique(m_temps)
    # u_ens = np.array([
    #     np.mean(m_ens[np.argwhere(m_temps == temp).flatten()]) for temp in u_temps]
    # )
    #
    # num_grad = np.diff(u_ens) / np.diff(u_temps)
    #
    # cp_func = np.poly1d(np.polyfit(u_temps[:-1], num_grad, 1))
    #
    # fig = go.Figure()
    # fig.add_scattergl(x=m_temps, y=m_ens, mode='markers')
    # fig.add_scattergl(x=u_temps, y=u_ens)
    # fig.add_scattergl(x=u_temps, y=num_grad, mode='markers')
    # fig.add_scattergl(x=u_temps, y=cp_func(u_temps))
    # fig.show(renderer='browser')

    'fit analytic gradient function'
    fit = np.polyfit(m_temps, m_ens, 2)
    melt_fit = np.poly1d(fit)

    fig = go.Figure()

    fig.add_scattergl(
        x=melt_enthalpies['Temperature'], y=melt_enthalpies['Enthalpies'],
        mode='markers', showlegend=False,
        marker_color=colors[-1],

    )

    fig.add_scattergl(
        x=melt_enthalpies['Temperature'], y=melt_fit(melt_enthalpies['Temperature']),
        mode='lines', showlegend=True, name='Melt ' + str(fit),
        marker_color=colors[-1],

    )

    for ident, (key, subdict) in enumerate(enthalpies_dict.items()):
        T = subdict['Temperature']
        H = subdict['Enthalpies']
        melt_T = melts_dict[key]

        solid_inds = np.argwhere((T < melt_T)).flatten()

        fig.add_scattergl(
            x=T[solid_inds], y=H[solid_inds],
            mode='markers', showlegend=False,
            marker_color=colors[ident],
        )

        fit2 = np.polyfit(T[solid_inds], H[solid_inds], 2)
        solid_fit = np.poly1d(fit2)
        xspace = np.linspace(T[solid_inds].min(), T[solid_inds].max(), 101)

        fig.add_scattergl(
            x=xspace, y=solid_fit(xspace),
            mode='lines', showlegend=True, name=key + ' ' + str(fit2),
            marker_color=colors[ident],

        )

        cp_function = solid_fit.deriv()
        latents_dict[key] = melt_fit(melt_T) - solid_fit(melt_T)
        cp_coefficients_dict[key] = cp_function.coefficients
        cp_at_T[key] = cp_function(melt_T)

    print(latents_dict)
    print(cp_coefficients_dict)
    print(cp_at_T)
    fig.update_xaxes(title='Temperature /K')
    fig.update_yaxes(title='Mean Enthalpy kJ/mol')
    fig.show(renderer='browser')

    fig = go.Figure()
    from e2emolmats.analysis.free_energy_calc import experimental_polymorphs_dict
    latents_dict2 = {key.replace('acridine/', ''): latents_dict[key] for key in latents_dict.keys()}
    exp_melts_dict = {key.replace(' ', ''): experimental_polymorphs_dict[key]['H_fus'] for key in
                      experimental_polymorphs_dict.keys()}
    fig.add_trace(
        go.Bar(x=list(latents_dict2.keys()), y=list(latents_dict2.values()), name='Latent Heat Estimate'))
    fig.add_trace(
        go.Bar(x=list(exp_melts_dict.keys()), y=list(exp_melts_dict.values()), name='Latent Heat (exp)'))
    fig.add_trace(go.Bar(x=['Reference 383K', 'Reference2 383K'], y=[20.682, 18.58], name='Reference'))
    fig.update_layout(yaxis_title='Latent Heat of Fusion (kJ/mol)')
    fig.show(renderer='browser')

    return fig


def extract_df_enthalpies(good_df, n_a):
    mean_enthalpy = np.zeros(len(good_df))
    temps = np.zeros_like(mean_enthalpy)

    for run_ind, row in good_df.iterrows():
        if 'E_tot' in row.keys():
            if np.sum(np.isnan(row['E_tot'])) == 0:
                en_key = 'E_tot'
            else:
                en_key = 'TotEng'
        else:
            en_key = 'TotEng'

        temps[run_ind] = row['temperature']
        energy = row[en_key] * 4.184 / row['num_molecules']  # / \
        pressure = 1  # good_df.iloc[run_ind]['Press']  # atmospheres
        volume = row['Volume']
        # atm*cubic angstrom = 1.01325e-28 kJ, normed per-mol
        # atm*cubic meter = 101.325 kJ = 101325 J
        volume_m3 = volume / 1e30
        PV_energy2 = pressure * volume_m3 * 101.325  # kilojoules per simulation box
        PV_energy = PV_energy2 / row['num_molecules'] * n_a  # kilojoules per mol
        #PV_energy = pressure * volume * n_a / row['num_molecules'] * (1.01325 * 10 ** -25) / 1000
        run_enthalpy = energy[-len(PV_energy):] + PV_energy

        equil_time = row['run_config']['equil_time']
        run_time = row['run_config']['run_time']

        crystal_reference_time = equil_time
        crystal_time_index = np.argmin(np.abs(row['time step'] - crystal_reference_time))

        melt_reference_time = 2 * equil_time
        melt_time_index = np.argmin(np.abs(row['time step'] - melt_reference_time))

        sampling_start_time = 5 * equil_time
        sampling_start_index = np.argmin(np.abs(row['time step'] - sampling_start_time))

        sampling_end_time = 5 * equil_time + run_time
        sampling_end_index = np.argmin(np.abs(row['time step'] - sampling_end_time))

        if row['Melt Succeeded']:
            mean_enthalpy[run_ind] = np.mean(run_enthalpy[sampling_start_index:])
        else:
            mean_enthalpy[run_ind] = np.mean(run_enthalpy[len(run_enthalpy) // 2:])

    return temps, mean_enthalpy


def relabel_defects(combined_df):
    """split 2,7-dihydroxynaphthalene sym and anti into two different defect types
    also, identify pure samples and relabel as such"""

    dhn_inds = np.argwhere(combined_df['defect_type'] == '2,7-dihydroxynaphthalene').flatten()
    defects_list = []
    for ind in range(len(combined_df)):
        if combined_df.iloc[ind]['defect_rate'] == 0:
            defects_list.append('Pure')
        elif ind in dhn_inds:
            inverted_defect = combined_df.iloc[ind]['run_config']['invert_defects']
            if inverted_defect:
                defects_list.append('2,7-dihydroxynaphthalene_anti')
            else:
                defects_list.append('2,7-dihydroxynaphthalene_sym')
        else:
            defects_list.append(combined_df.iloc[ind]['defect_type'])

    combined_df['defect_type'] = defects_list

    return combined_df

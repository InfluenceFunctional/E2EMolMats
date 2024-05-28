import io
import os

import numpy as np
from _plotly_utils.colors import n_colors
from scipy.optimize import minimize_scalar
from scipy.stats import linregress
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from scipy.ndimage import gaussian_filter1d
import plotly.express as px


def make_thermo_fig(traj_thermo_keys, thermo_results_dict, run_config):
    thermo_telemetry_fig = make_subplots(rows=2, cols=5, subplot_titles=traj_thermo_keys)
    ind = 0
    for i, key in enumerate(thermo_results_dict.keys()):
        if key in traj_thermo_keys:
            row = ind // 5 + 1
            col = ind % 5 + 1
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
    good_keys = ['box_type', 'min_inter_cluster_distance', 'pressure_direction', 'run_name', 'damping',
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
    return thermo_telemetry_fig, thermo_results_dict['thermo_trajectory'].shape[1]


def process_thermo_data():
    results_dict = {'time step': [],
                    'temp': [],
                    'E_pair': [],
                    'E_mol': [],
                    'E_tot': [],
                    'PotEng': [],
                    'Press': [],
                    'Volume': [],
                    }
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
                skip = False
                # print(ind)
        else:
            if "Loop" in line:
                skip = True

            if not skip:
                split_line = line.split(' ')
                entries = [float(entry) for entry in split_line if entry != '']
                for ind2, key in enumerate(results_dict.keys()):
                    if ind2 < len(entries):
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
        return results_dict, 'Thermo analysis succeeded'

    else:
        return results_dict, 'Missing thermo analysis!'


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
    num_polymorphs = len(polymorphs)
    colors = px.colors.qualitative.G10
    seen_polymorph = {polymorph: False for polymorph in polymorphs}
    min_temp = np.amin([df.iloc[ind]['run_config']['temperature'] for ind in range(len(df))])
    max_temp = np.amax([df.iloc[ind]['run_config']['temperature'] for ind in range(len(df))])
    temprange = np.linspace(min_temp, max_temp, 1000)
    melt_temps = {}
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=['Normed Intermolecular Energy', 'Intermolecular Energy Slope'])
    for polymorph in polymorphs:
        good_inds = np.argwhere(df['polymorph_name'] == polymorph).flatten()

        temperatures = np.asarray([elem['temperature'] for elem in df.iloc[good_inds]['run_config']])
        melt_magnitudes = np.asarray(df.iloc[good_inds]['melt_magnitude']).flatten()
        melt_slopes = np.asarray(df.iloc[good_inds]['melt_slope']).flatten()

        temps = np.unique(temperatures)
        mag_at_t = np.array(
            [np.mean([melt_magnitudes[i] for i in range(len(melt_magnitudes)) if temperatures[i] == temp]) for temp
             in
             temps])
        slope_at_t = np.array(
            [np.mean([melt_slopes[i] for i in range(len(melt_slopes)) if temperatures[i] == temp]) for temp in
             temps])

        mag_spline = np.interp(temprange, temps, np.maximum.accumulate(mag_at_t))
        slope_spline = np.interp(temprange, temps, np.maximum.accumulate(slope_at_t))

        melt_T = temprange[np.argmin(np.abs(slope_spline))]
        melt_temps[polymorph] = melt_T
        fig.add_scattergl(x=temperatures,
                          y=melt_magnitudes,
                          mode='markers',
                          name=polymorph,
                          legendgroup=polymorph,
                          marker_size=7,
                          opacity=0.5,
                          showlegend=False,
                          marker_color=colors[polymorphs.index(polymorph)],
                          row=1, col=1
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
                          row=1, col=2
                          )

        fig.add_scattergl(x=temps,
                          y=mag_at_t,
                          mode='markers',
                          name=polymorph,
                          legendgroup=polymorph,
                          marker_size=15,
                          showlegend=True if not seen_polymorph[polymorph] else False,
                          marker_color=colors[polymorphs.index(polymorph)],
                          row=1, col=1
                          )

        fig.add_scattergl(x=temps,
                          y=slope_at_t,
                          mode='markers',
                          name=polymorph,
                          legendgroup=polymorph,
                          marker_size=15,
                          showlegend=False,
                          marker_color=colors[polymorphs.index(polymorph)],
                          row=1, col=2
                          )
        #
        # fig.add_scattergl(x=temprange,
        #                   y=mag_spline,
        #                   name=polymorph,
        #                   legendgroup=polymorph,
        #                   showlegend=False,
        #                   marker_color=colors[polymorphs.index(polymorph)],
        #                   row=1, col=1
        #                   )
        #
        # fig.add_scattergl(x=temprange,
        #                   y=slope_spline,
        #                   name=polymorph,
        #                   legendgroup=polymorph,
        #                   showlegend=False,
        #                   marker_color=colors[polymorphs.index(polymorph)],
        #                   row=1, col=2
        #                   )
    fig.update_xaxes(title='Temperature (K)')
    if show_fig:
        fig.show(renderer='browser')
    return fig, melt_temps


def plot_melt_points(melt_estimate_dict, true_melts_dict, show_fig=True):
    """
    use all data to estimate the melt point for each polymorph
    """
    melt_temps = np.asarray(list(melt_estimate_dict.values()))
    true_melt_temps = np.asarray(list(true_melts_dict.values()))
    fig2 = go.Figure()
    fig2.add_trace(go.Bar(x=list(melt_estimate_dict.keys()), y=list(melt_estimate_dict.values()), name='Estimate'))
    fig2.add_trace(go.Bar(x=list(true_melts_dict.keys()), y=list(true_melts_dict.values()), name='Reference'))
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


def runs_summary_table(runs_dict):
    runs = np.asarray(list(runs_dict.keys())).astype(int)
    runs = np.sort(runs)
    codes = [runs_dict[str(key)][0] for key in runs]
    defects = [runs_dict[str(key)][1]['defect_type'] for key in runs]
    defect_rates = [runs_dict[str(key)][1]['defect_rate'] for key in runs]
    sizes = [runs_dict[str(key)][1]['max_sphere_radius'] for key in runs]
    polymorph = [runs_dict[str(key)][1]['structure_identifier'].split('/')[-1] for key in runs]

    fig = go.Figure(
        data=go.Table(header=dict(
            values=['run_num', 'run_status', 'defect_type', 'defect_rate', 'cluster_radius', 'polymorph']),
            cells=dict(values=[runs, codes, defects, defect_rates, sizes, polymorph])))
    return fig
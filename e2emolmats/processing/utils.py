import io
import os

import numpy as np
from plotly import graph_objects as go
from plotly.subplots import make_subplots
from scipy.ndimage import gaussian_filter1d
from scipy.stats import linregress

from e2emolmats.reporting.temp_analysis import extract_local_profile
from e2emolmats.reporting.utils import extract_trajectory_properties, read_lammps_thermo_traj, read_lammps_com_traj, \
    com_mobility, read_lammps_pe_traj


def mode_settings(MODE,
                  acridine_cluster_paths,
                  acridine_melt_paths,
                  acridine_scan_paths,
                  acridine_latent_paths,
                  acridine_cp_paths,
                  acridine_cp2_paths,
                  acridine_lattice_energy_paths,
                  ):
    compute_melt_temps = False
    nanocluster_analysis = False
    latents_analysis = False
    cp_analysis = False
    cp2_analysis = False
    lattice_energy_analysis = False
    melt_scan_analysis = False

    if MODE == 'acridine_cluster':
        battery_paths = acridine_cluster_paths
        nanocluster_analysis = True

    elif MODE == 'acridine_melt':
        battery_paths = acridine_melt_paths
        compute_melt_temps = True

    elif MODE == 'acridine_scan':
        battery_paths = acridine_scan_paths
        melt_scan_analysis = True

    elif MODE == 'acridine_latent':
        battery_paths = acridine_latent_paths
        latents_analysis = True

    elif MODE == 'acridine_cp':
        battery_paths = acridine_cp_paths
        cp_analysis = True

    elif MODE == 'acridine_cp2':
        battery_paths = acridine_cp2_paths
        cp2_analysis = True

    elif MODE == 'acridine_lattice_energy':
        battery_paths = acridine_lattice_energy_paths
        lattice_energy_analysis = True

    else:
        assert False, "Unrecognized mode !!"

    return (battery_paths,
            melt_scan_analysis,
            nanocluster_analysis,
            compute_melt_temps,
            latents_analysis,
            cp_analysis,
            cp2_analysis,
            lattice_energy_analysis
            )


def make_thermo_figs(thermo_results_dict, run_config):
    keys_to_use = []
    for key in thermo_results_dict.keys():
        if isinstance(thermo_results_dict['Temp'], np.ndarray):
            if thermo_results_dict[key].ndim == 1:
                keys_to_use += [key]

    good_keys = ['box_type', 'min_inter_cluster_distance', 'run_name', 'damping',
                 'defect_rate', 'defect_type', 'gap_rate', 'max_sphere_radius', 'min_lattice_length',
                 'scramble_rate', 'seed', 'structure_identifier', 'temperature']

    #thermo_telemetry_fig = make_subplots(rows=int(np.ceil(n_keys / cols)), cols=cols, subplot_titles=keys_to_use)
    thermo_figs_dict = {}
    ind = 0
    for i, key in enumerate(keys_to_use):
        fig = go.Figure()
        fig.add_trace(
            go.Scattergl(x=thermo_results_dict['time step'] / 1e6,
                         y=thermo_results_dict[key],
                         name=key, showlegend=False),
        )
        fig.add_trace(
            go.Scattergl(x=thermo_results_dict['time step'] / 1e6,
                         y=gaussian_filter1d(thermo_results_dict[key], sigma=2),
                         name=key, showlegend=False),
        )
        fig.update_layout(xaxis_title="Time (ns)", yaxis_title=keys_to_use[i])
        fig.update_layout(title=str({key: run_config[key] for key in good_keys[0:5]})
                                + '<br>' + str({key: run_config[key] for key in good_keys[5:10]})
                                + '<br>' + str({key: run_config[key] for key in good_keys[10:]}),
                          title_automargin=True,
                          title_pad_t=20,
                          margin_t=180)

        thermo_figs_dict[key] = fig
        ind += 1

    profile_keys = []
    for key in thermo_results_dict.keys():
        if thermo_results_dict[key].ndim == 2 and 'profile' in key:
            profile_keys.append(key)

    for ind, key in enumerate(profile_keys):
        fig = go.Figure()
        fig.add_heatmap(z=thermo_results_dict[key][thermo_results_dict['sampling_start_index']:], colorbar_title=key)
        fig.update_layout(title=str({key: run_config[key] for key in good_keys[0:5]})
                                + '<br>' + str({key: run_config[key] for key in good_keys[5:10]})
                                + '<br>' + str({key: run_config[key] for key in good_keys[10:]}),
                          title_automargin=True,
                          title_pad_t=20,
                          margin_t=180)
        thermo_figs_dict[key] = fig

    # thermo_telemetry_fig.update_layout(title=f"{run_config['structure_identifier']}, "
    #                                          f"Gap Rate {run_config['gap_rate']:.2f}, "
    #                                          f"System Size {run_config['min_lattice_length']} "
    #                                          f"Cubic Angstrom or {thermo_results_dict['thermo_trajectory'].shape[1]} Molecules"
    #                                          f"Sampling temperature {run_config['temperature']}K"
    #                                    )
    return thermo_figs_dict


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
        if 'Target temperature for fix nvt cannot be 0.0' in text \
                or 'Target temperature for fix npt cannot be 0.0' in text:  # minimization only
            results_dict['total time'] = 0
        else:
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

    results_dict = extract_trajectory_properties(hit_minimization, lines, results_dict, skip)

    delta_t = np.mean(np.diff(results_dict['Step'][-20:]))
    first_delta_t_index = np.argmin(np.abs(delta_t - np.array(results_dict['Step'])))
    results_dict['time step'] = np.array(
        results_dict['Step'][max(first_delta_t_index - 1, 0):])  # correct for minimizations

    # label changes between simulation fixes
    equil_time = run_config['equil_time']

    if (run_config['run_name'] == 'acridine_melt_interface5' or
            run_config['run_name'] == 'acridine_interface_scan2' or
            run_config['run_name'] == 'acridine_interface_scan3'):
        melt_time = 100000
    else:
        melt_time = equil_time
    sampling_time = run_config['run_time']

    init_equil_time = equil_time
    init_equil_time_index = np.argmin(np.abs(results_dict['time step'] - init_equil_time))
    melt_time = equil_time + melt_time
    melt_time_index = np.argmin(np.abs(results_dict['time step'] - melt_time))
    sampling_start_time = 4 * equil_time + melt_time
    sampling_start_index = np.argmin(np.abs(results_dict['time step'] - sampling_start_time))
    sampling_end_time = sampling_time + sampling_start_time
    sampling_end_index = np.argmin(np.abs(results_dict['time step'] - sampling_end_time))

    fix_label = np.zeros_like(results_dict['time step'])
    for ind in range(len(fix_label)):
        if ind < init_equil_time_index:
            pass  # init equilibration
        elif init_equil_time_index <= ind < melt_time_index:
            fix_label[ind] = 1  # melting
        elif melt_time_index <= ind < sampling_start_index:
            fix_label[ind] = 2  # cooling
        elif sampling_start_index <= ind < sampling_end_index:
            fix_label[ind] = 3  # sampling
        else:
            fix_label[ind] = 4  # post-sampling sampling

    results_dict['fix_label'] = fix_label
    results_dict.update({
        'equil_time_index': init_equil_time_index,
        'melt_time_index': melt_time_index,
        'sampling_start_index': sampling_start_index,
        'sampling_end_index': sampling_end_index
    })

    for key in results_dict.keys():
        results_dict[key] = np.asarray(results_dict[key])

    if not skip_molwise_thermo:
        if os.path.exists('tmp.out'):  # molecule-wise temperature analysis
            frames = read_lammps_thermo_traj('tmp.out')
            results_dict['thermo_trajectory'] = np.asarray(list(frames.values()))

        if os.path.exists('com.out'):  # molecule-wise temperature analysis
            frames = read_lammps_com_traj('com.out')
            results_dict['com_trajectory'] = np.asarray(list(frames.values()))
            results_dict['com_mobility'] = com_mobility(results_dict['com_trajectory'],
                                                        trailing_length=5)

        if os.path.exists('pe.out'):  # molecule-wise temperature analysis
            frames = read_lammps_pe_traj('pe.out')
            results_dict['pe_trajectory'] = np.asarray(list(frames.values()))

        if os.path.exists('pe.out') and os.path.exists('com.out') and os.path.exists('tmp.out'):
            # if we have all thermo outputs, process several extra outputs

            """
            for property
                pe, mobility, T1-T3
                for location
                    all bulk, all melt, interface bulk, interface melt, deep bulk, deep melt
                        log 1d feature
            """
            props_list = ['Mol T', 'KE Com', 'Intra T', 'Mobility', 'PE']
            all_props_traj = np.concatenate([
                results_dict['thermo_trajectory'],
                results_dict['com_mobility'][..., None],
                results_dict['pe_trajectory']
            ],
                axis=2)

            'Full bulk and crystal properties'
            types_list = ['Bulk', 'Melt']
            crystal_inds = np.arange(run_config['melt_indices'].crystal_start_ind,
                                     run_config['melt_indices'].crystal_end_ind)
            melt_inds = np.arange(run_config['melt_indices'].melt_start_ind,
                                  run_config['melt_indices'].melt_end_ind)

            inds_list = [crystal_inds, melt_inds]
            for ind1, type in enumerate(types_list):
                for ind2, prop_name in enumerate(props_list):
                    key = f"{type} {prop_name}"
                    prop = all_props_traj[:, inds_list[ind1], ind2].mean(-1)
                    results_dict[key] = prop

            'Local property profiles'
            nbinsx, nbinsy = 30, len(results_dict['com_trajectory']) // 2
            interface1_ind = nbinsx // 2
            interface2_ind = 0
            mid_bulk_ind = nbinsx // 4
            mid_melt_ind = nbinsx // (3 / 2)

            location_inds = {
                'interface1': interface1_ind,
                'interface2': interface2_ind,
                'mid_bulk': mid_bulk_ind,
                'mid_melt': mid_melt_ind,
            }

            for ind2, prop_name in enumerate(props_list):
                prop = all_props_traj[..., ind2]
                for compute_anomaly in [True, False]:
                    if compute_anomaly:
                        key = f"{prop_name} anomaly profile"
                    else:
                        key = f"{prop_name} profile"

                    profile, bins_x, bins_y = extract_local_profile(
                        run_config,
                        results_dict['time step'],
                        prop,
                        results_dict['com_trajectory'],
                        nbinsx, nbinsy,
                        compute_anomaly=compute_anomaly
                    )
                    results_dict[key] = profile
                    for l_ind, (location_name, location) in enumerate(location_inds.items()):
                        key += f' {location_name}'
                        results_dict[key] = profile[:, int(location)]

            return results_dict, 'Thermo analysis succeeded'
        else:

            return results_dict, 'Missing thermo trajectory!'
    else:
        return results_dict, 'Thermo analysis succeeded'


def get_melt_progress(results_df):
    assert False, "Rewrite this using CoM mobility"
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

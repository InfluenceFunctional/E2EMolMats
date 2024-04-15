import os
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from e2emolmats.reporting.utils import process_thermo_data, make_thermo_fig, multi_ramp_fig, make_gap_vs_tm_fig
from e2emolmats.common.utils import dict2namespace
import pandas as pd
import wandb
import glob

battery_paths = [
    r'D:\crystal_datasets\acridine_melt_interface1/',
]

combined_df = pd.DataFrame()

for battery_path in battery_paths:
    params = {
        'battery_path': battery_path,
        'machine': 'local',  # or 'cluster'  ### doesn't do anything
        'show_figs': False,
        'make_run_wise_figs': True,
    }
    config = dict2namespace(params)

    wandb.init(config=params, project="E2EMolMats",
               entity="mkilgour", tags=[config.battery_path],
               settings=wandb.Settings(code_dir="../../analysis"))

    wandb.run.name = config.battery_path
    wandb.run.save()

    os.chdir(config.battery_path)
    battery_full_path = os.getcwd()

    if os.path.exists('results_df'):
        results_df = pd.read_pickle('results_df')
    else:
        results_df = pd.DataFrame(columns=["run_num",
                                           "pressure",
                                           "E_pair",
                                           "E_mol",
                                           "E_tot",
                                           ])

    traj_thermo_keys = ['temp', 'E_pair', 'E_mol', 'E_tot', 'Press', 'Volume', 'molwise_mean_temp',
                        'molwise_mean_kecom', 'molwise_mean_internal']

    dirs = os.listdir()
    dirs += glob.glob('*/*')
    for run_dir in dirs:  # loop over run directories in the battery
        os.chdir(config.battery_path)
        # (run_dir not in results_df["run_num"].values) and \

        if (run_dir != 'md_data') and \
                (run_dir not in results_df["run_num"].values) and \
                ('results_df' not in run_dir) and \
                ('png' not in run_dir) and \
                ('log' not in run_dir) and \
                ('wandb' not in run_dir):

            try:
                os.chdir(run_dir)
                # do the analysis
                print(run_dir)
                run_config = np.load('run_config.npy', allow_pickle=True).item()
            except:
                continue

            if config.make_run_wise_figs:
                '''thermodynamic data'''
                thermo_results_dict = process_thermo_data()
                if 'thermo_trajectory' not in thermo_results_dict.keys():
                    continue
                thermo_telemetry_fig = make_thermo_fig(traj_thermo_keys, thermo_results_dict, run_config)

                if config.show_figs:
                    thermo_telemetry_fig.show(renderer="browser")
                wandb.log({'Thermo Data': thermo_telemetry_fig})
                # thermo_telemetry_fig.write_image('Thermo Data.png')

                '''save results'''
                new_row = {"run_num": run_dir,
                           "gap_rate": [run_config['gap_rate']],
                           'min_lattice_length': [run_config['min_lattice_length']],
                           'num_molecules': [thermo_results_dict['thermo_trajectory'].shape[1]],
                           'run_config': [run_config],
                           }
                for key in traj_thermo_keys:
                    new_row.update({key: [thermo_results_dict[key]]})
                results_df = pd.concat([results_df, pd.DataFrame.from_dict(new_row)])
                results_df.to_pickle(battery_full_path + '/results_df')

    results_df.reset_index(drop=True, inplace=True)

    # multi-run analysis
    '''combined temperature ramps'''
    temp_vs_pe_fig = multi_ramp_fig(results_df)

    if config.show_figs:
        temp_vs_pe_fig.show(renderer='browser')
    wandb.log({"Temperature Ramp": temp_vs_pe_fig})

    '''defect rate vs melt point, per-polymorph'''
    # gap_vs_tm_fig, results_df, melt_fit_figs = make_gap_vs_tm_fig(results_df)
    #
    # if config.show_figs:
    #     gap_vs_tm_fig.show(renderer='browser')
    #
    # wandb.log({"Gap rate vs T Melt": gap_vs_tm_fig})
    if len(combined_df) > 0:
        combined_df = pd.concat([combined_df, results_df])
    else:
        combined_df = results_df

    wandb.finish()

'''
combination of several series

plot melting point vs run length, system size
'''
combined_df.reset_index(drop=True, inplace=True)

runtimes = [df['run_config']['run_time'] / 1e6 for _, df in combined_df.iterrows()]

polymorphs = list(np.unique([conf['structure_identifier'].split('/')[-1] for conf in combined_df['run_config']]))
num_polymorphs = len(polymorphs)

import plotly.express as px

colors = px.colors.qualitative.G10
seen_polymorph = {polymorph: False for polymorph in polymorphs}

fig = go.Figure()
for polymorph in polymorphs:
    good_inds = np.argwhere((combined_df['polymorph_name'] == polymorph) * (combined_df['step_size'] > 0.25)).flatten()
    fig.add_scattergl(x=np.asarray(runtimes)[good_inds],
                      y=np.asarray(combined_df.iloc[good_inds]['polymorph_melt_temp']).flatten(),
                      mode='markers',
                      name=polymorph,
                      legendgroup=polymorph,
                      # marker_size=15,
                      showlegend=True if not seen_polymorph[polymorph] else False,
                      marker_color=colors[polymorphs.index(polymorph)],
                      )

fig.update_layout(xaxis_title='Run Time (ns)', yaxis_title='Tm (K)')
fig.show(renderer='browser')
fig.write_image('../../melt_vs_runtime.png')
np.save('../../combined_analysis_dict', combined_df)
'''big combined multivariate figure'''

fig = make_subplots(rows=1, cols=3)
for polymorph in polymorphs:
    good_inds = np.argwhere((combined_df['polymorph_name'] == polymorph) * (combined_df['step_size'] > 0.25)).flatten()
    y = np.asarray(combined_df.iloc[good_inds]['melt_temp']).flatten()
    x1 = np.asarray(runtimes)[good_inds]
    x1u = np.unique(x1)
    y1_means = [np.mean(y[np.argwhere(x1 == unique).flatten()]) for unique in x1u]
    x2 = np.asarray(combined_df.iloc[good_inds]['min_lattice_length'])
    x2u = np.unique(x2)
    y2_means = [np.mean(y[np.argwhere(x2 == unique).flatten()]) for unique in x2u]
    x3 = np.asarray(combined_df.iloc[good_inds]['gap_rate'])
    x3u = np.unique(x3)
    y3_means = [np.mean(y[np.argwhere(x3 == unique).flatten()]) for unique in x3u]
    fig.add_scattergl(x=x1,
                      y=y,
                      mode='markers',
                      name=polymorph,
                      legendgroup=polymorph,
                      showlegend=True if not seen_polymorph[polymorph] else False,
                      marker_color=colors[polymorphs.index(polymorph)],
                      row=1, col=1
                      )

    fig.add_scattergl(x=x2,
                      y=np.asarray(combined_df.iloc[good_inds]['melt_temp']).flatten(),
                      mode='markers',
                      name=polymorph,
                      legendgroup=polymorph,
                      showlegend=False,  #True if not seen_polymorph[polymorph] else False,
                      marker_color=colors[polymorphs.index(polymorph)],
                      row=1, col=2
                      )
    fig.add_scattergl(x=x3,
                      y=np.asarray(combined_df.iloc[good_inds]['melt_temp']).flatten(),
                      mode='markers',
                      name=polymorph,
                      legendgroup=polymorph,
                      showlegend=False,  #True if not seen_polymorph[polymorph] else False,
                      marker_color=colors[polymorphs.index(polymorph)],
                      row=1, col=3
                      )
    fig.add_scattergl(x=x1u,
                      y=y1_means,
                      mode='markers',
                      name=polymorph,
                      legendgroup=polymorph,
                      showlegend=False,
                      marker_size=15,
                      marker_color=colors[polymorphs.index(polymorph)],
                      row=1, col=1
                      )
    fig.add_scattergl(x=x2u,
                      y=y2_means,
                      mode='markers',
                      name=polymorph,
                      legendgroup=polymorph,
                      showlegend=False,
                      marker_size=15,
                      marker_color=colors[polymorphs.index(polymorph)],
                      row=1, col=2
                      )
    fig.add_scattergl(x=x3u,
                      y=y3_means,
                      mode='markers',
                      name=polymorph,
                      legendgroup=polymorph,
                      showlegend=False,
                      marker_size=15,
                      marker_color=colors[polymorphs.index(polymorph)],
                      row=1, col=3
                      )

fig.update_yaxes(title='Melt Temperature (K)', row=1, col=1)
fig.update_xaxes(title='Runtime (ns)', row=1, col=1)
fig.update_xaxes(title='Supercell Edge Length (A)', row=1, col=2)
fig.update_xaxes(title='Gap Fraction', row=1, col=3)

fig.show(renderer='browser')
fig.write_image('../../melt_vs_runtime.png')

fig = go.Figure()
for polymorph in ["Form2"]:
    good_inds = np.argwhere((combined_df['polymorph_name'] == polymorph) * (combined_df['step_size'] > 0.25)).flatten()
    y = np.asarray(combined_df.iloc[good_inds]['melt_temp']).flatten()

    x1 = np.asarray(runtimes)[good_inds]
    x1u = np.unique(x1)

    x3 = np.asarray(combined_df.iloc[good_inds]['gap_rate'])
    x3u = np.unique(x3)

    y_means = np.asarray(
        [np.mean(y[np.argwhere((x1 == unique1) * (x3 == unique2)).flatten()]) for unique1 in x1u for unique2 in x3u])
    x1ul = np.asarray([unique1 for unique1 in x1u for unique2 in x3u])
    x3ul = np.asarray([unique2 for unique1 in x1u for unique2 in x3u])

    good_inds = np.argwhere(np.isfinite(y_means)).flatten().astype(int)
    fig.add_scattergl(x=x3ul[good_inds], y=x1ul[good_inds], marker_color=y_means[good_inds], mode='markers', opacity=1,
                      marker_colorbar=dict(title="Melt Temperature (K)"))

fig.update_layout(xaxis_title='Gap Rate', yaxis_title='Runtime (ns)')

fig.show(renderer='browser')
fig.write_image('../../gap_vs_runtime_vs_melt.png')
aa = 1

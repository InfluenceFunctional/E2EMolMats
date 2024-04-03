import os
import numpy as np
import plotly.graph_objects as go

from e2emolmats.reporting.utils import process_thermo_data, make_thermo_fig, multi_ramp_fig, make_gap_vs_tm_fig
from e2emolmats.common.utils import dict2namespace
import pandas as pd
import wandb
import glob

battery_paths = [
    # r'D:\crystal_datasets\acridine_melt_series2_1/',
    # r'D:\crystal_datasets\acridine_melt_series2_2/',
    # r'D:\crystal_datasets\acridine_melt_series2_3/',
    #r'D:\crystal_datasets\acridine_melt_series2_4/',
    #r'D:\crystal_datasets\acridine_melt_series2_5/',
    r'D:\crystal_datasets\acridine_melt_series2_6/'
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
    gap_vs_tm_fig, results_df = make_gap_vs_tm_fig(results_df)

    if config.show_figs:
        gap_vs_tm_fig.show(renderer='browser')

    wandb.log({"Gap rate vs T Melt": gap_vs_tm_fig})
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
    good_inds = np.argwhere(combined_df['polymorph_name'] == polymorph).flatten()
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

aa = 1

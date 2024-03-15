import os
import numpy as np
from reporting.utils import process_thermo_data
from utils import dict2namespace
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import wandb

params = {
    'battery_path': r'D:\crystal_datasets\acridine_pure_bulk_melt2/',
    'machine': 'local',  # or 'cluster'  ### doesn't do anything
    'show_figs': False,
    'make_run_wise_figs': True,
}
config = dict2namespace(params)

wandb.init(config=params, project="nicotinamide_clusters",
           entity="mkilgour", tags=["acridine_pure_bulk_melt2"],
           settings=wandb.Settings(code_dir="."))

wandb.run.name = config.battery_path
wandb.run.save()

os.chdir(config.battery_path)

if os.path.exists('results_df'):
    results_df = pd.read_pickle('results_df')
else:
    results_df = pd.DataFrame(columns=["run_num",
                                       "pressure",
                                       "E_pair",
                                       "E_mol",
                                       "E_tot",
                                       ])
#
# dirs = os.listdir()
# for run_dir in dirs:  # loop over run directories in the battery
#     os.chdir(config.battery_path)
#     # (run_dir not in results_df["run_num"].values) and \
#
#     if (run_dir != 'common') and \
#             ('results_df' not in run_dir) and \
#             ('png' not in run_dir):
#         os.chdir(run_dir)
#         # do the analysis
#         print(run_dir)
#
#         run_config = np.load('run_config.npy', allow_pickle=True).item()
#
#         if config.make_run_wise_figs:
#             '''thermodynamic data'''
#             traj_thermo_keys = ['temp', 'E_pair', 'E_mol', 'E_tot', 'Press', 'molwise_mean_temp',
#                                 'molwise_mean_kecom', 'molwise_mean_internal']
#             thermo_results_dict = process_thermo_data()
#             fig = make_subplots(rows=2, cols=4, subplot_titles=traj_thermo_keys)
#             ind = 0
#             for i, key in enumerate(thermo_results_dict.keys()):
#                 if key in traj_thermo_keys:
#                     row = ind // 4 + 1
#                     col = ind % 4 + 1
#                     fig.add_trace(
#                         go.Scattergl(x=thermo_results_dict['time step'] / 1e6,
#                                      y=thermo_results_dict[key],
#                                      name=key, showlegend=False),
#                         row=row, col=col
#                     )
#                     ind += 1
#
#             fig.update_xaxes(title_text="Time (ns)")
#             fig.update_layout(title=f"{run_config['structure_identifier']}, "
#                                     f"Gap Rate {run_config['gap_rate']:.2f}, "
#                                     f"System Size {run_config['min_lattice_length']} Cubic Angstrom or {thermo_results_dict['thermo_trajectory'].shape[1]} Molecules"
#                               )
#
#             if config.show_figs:
#                 fig.show(renderer="browser")
#             wandb.log({'Thermo Data': fig})
#             fig.write_image('Thermo Data.png')
#
#             '''save results'''
#             new_row = {"run_num": run_dir,
#                        "gap_rate": [run_config['gap_rate']],
#                        'min_lattice_length': [run_config['min_lattice_length']],
#                        'num_molecules': [thermo_results_dict['thermo_trajectory'].shape[1]],
#                        }
#             for key in traj_thermo_keys:
#                 new_row.update({key: [thermo_results_dict[key]]})
#             results_df = pd.concat([results_df, pd.DataFrame.from_dict(new_row)])
#             results_df.to_pickle('../results_df')

# multi-run analysis
results_df.reset_index(drop=True, inplace=True)
from plotly.colors import n_colors

colors = n_colors('rgb(5,120,200)', 'rgb(250,50,5)', len(np.unique(results_df['gap_rate'])), colortype='rgb')
gaps = np.sort(np.unique(results_df['gap_rate']))
gap_dict = {temp: ind for ind, temp in enumerate(gaps)}

fig = go.Figure()
for ind in range(len(results_df)):
    gap = results_df['gap_rate'][ind]
    fig.add_scattergl(x=results_df['temp'][ind],
                      y=results_df['E_pair'][ind],
                      line_color=colors[gap_dict[gap]],
                      marker_color=gap,
                      mode='markers',
                      marker_size=5,
                      showlegend=False,
                      opacity=0.5,
                      # name=f"Gap rate = {gap:.2f}",
                      # showlegend=True
                      )
fig.update_layout(xaxis_title="Temperature (K)", yaxis_title="Pair Energy")
if config.show_figs:
    fig.show(renderer='browser')

wandb.log({"Temperature Ramp": fig})
aa = 1

import os
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats import linregress
import plotly.express as px

from e2emolmats.reporting.utils import process_thermo_data, make_thermo_fig, get_melt_progress
from e2emolmats.common.utils import dict2namespace
import pandas as pd
import wandb
import glob

if __name__ == '__main__':
    battery_paths = [
        r'D:\crystal_datasets\acridine_melt_interface14/',
        r'D:\crystal_datasets\acridine_melt_interface15/',
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

        traj_thermo_keys = ['temp', 'E_pair', 'E_mol', 'E_tot', 'PotEng',
                            'Press', 'Volume', 'molwise_mean_temp',
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

                    '''save results'''
                    new_row = {"run_num": run_dir,
                               "gap_rate": [run_config['gap_rate']],
                               'min_lattice_length': [run_config['min_lattice_length']],
                               'num_molecules': [thermo_results_dict['thermo_trajectory'].shape[1]],
                               'run_config': [run_config],
                               }
                    for key in traj_thermo_keys:
                        new_row.update({key: [thermo_results_dict[key]]})
                    new_row.update({'time step': [thermo_results_dict['time step']]})
                    results_df = pd.concat([results_df, pd.DataFrame.from_dict(new_row)])
                    results_df.to_pickle(battery_full_path + '/results_df')

        results_df.reset_index(drop=True, inplace=True)
        results_df['melt_slope'], results_df['melt_magnitude'] = get_melt_progress(results_df)
        results_df.to_pickle(battery_full_path + '/results_df')

        wandb.finish()
        if len(combined_df) > 0:
            combined_df = pd.concat([combined_df, results_df])
        else:
            combined_df = results_df

    polymorph_names = []
    seeds = []
    for _, row in combined_df.iterrows():
        polymorph_names.append(row['run_config']['structure_identifier'].split('/')[1])
        seeds.append(row['run_config']['seed'])
    combined_df['polymorph_name'] = polymorph_names
    combined_df['seed'] = seeds
    combined_df.reset_index(drop=True, inplace=True)

    seeds = list(np.unique(combined_df['seed']))
    polymorphs = list(np.unique([conf['structure_identifier'].split('/')[-1] for conf in combined_df['run_config']]))
    num_polymorphs = len(polymorphs)

    colors = px.colors.qualitative.G10
    seen_polymorph = {polymorph: False for polymorph in polymorphs}
    min_temp = np.amin([combined_df.iloc[ind]['run_config']['temperature'] for ind in range(len(combined_df))])
    max_temp = np.amax([combined_df.iloc[ind]['run_config']['temperature'] for ind in range(len(combined_df))])

    temprange = np.linspace(min_temp, max_temp, 1000)
    melt_temps = {}
    fig = make_subplots(rows=1, cols=2, subplot_titles=['Normed Intermolecular Energy', 'Intermolecular Energy Slope'])
    for polymorph in polymorphs:
        good_inds = np.argwhere(combined_df['polymorph_name'] == polymorph).flatten()

        temperatures = np.asarray([elem['temperature'] for elem in combined_df.iloc[good_inds]['run_config']])
        melt_magnitudes = np.asarray(combined_df.iloc[good_inds]['melt_magnitude']).flatten()
        melt_slopes = np.asarray(combined_df.iloc[good_inds]['melt_slope']).flatten()

        temps = np.unique(temperatures)
        mag_at_t = np.array([np.mean([melt_magnitudes[i] for i in range(len(melt_magnitudes)) if temperatures[i] == temp]) for temp in temps])
        slope_at_t = np.array([np.mean([melt_slopes[i] for i in range(len(melt_slopes)) if temperatures[i] == temp]) for temp in temps])

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
    fig.show(renderer='browser')

    '''
    use all data to estimate the melt point for each polymorph
    '''
    fig = go.Figure()
    fig.add_trace(go.Bar(x=list(melt_temps.keys()), y=list(melt_temps.values())))
    fig.show(renderer='browser')

    aa = 1

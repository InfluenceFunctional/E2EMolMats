import MDAnalysis as mda
import os
import warnings
from random import shuffle

warnings.filterwarnings("ignore", category=DeprecationWarning)  # ignore numpy error

from reporting.cluster_figs import (plot_thermodynamic_data, trajectory_rdf_analysis,
                                    cluster_molecule_alignment, process_thermo_data)
from utils import (dict2namespace, rewrite_trajectory, compute_Ip_molwise_alignment, process_dump)
import numpy as np

import pandas as pd
import plotly.io as pio
import tqdm

pio.renderers.default = 'browser'

# params = {
#     'reference_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\bulk_reference/',
#     'battery_path': r'D:\crystals_extra\defect_clusters_5_rerun/',
#     'machine': 'local',  # or 'cluster'  ### doesn't do anything
#     'show_figs': False,
#     'write_trajectory': False,
#     'make_run_wise_figs': True,
#     'do_rdf_analysis': True,
#     'do_alignment_analysis': True,
#     'results_df_path': 'results_df',
#     'reference_df_path': 'reference_df5',
#     'do_reference_analysis': False,
#     'do_sample_analysis': True,
#     'do_NN_analysis': True
# }

params = {
    'reference_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\bulk_reference/',
    'battery_path': r'/vast/mk8347/molecule_clusters/defect_clusters_6',
    'machine': 'cluster',  # or 'cluster'  ### doesn't do anything
    'show_figs': False,
    'write_trajectory': False,
    'make_run_wise_figs': True,
    'do_rdf_analysis': True,
    'do_alignment_analysis': True,
    'results_df_path': 'results_df',
    'reference_df_path': 'reference_df5',
    'do_reference_analysis': False,
    'do_sample_analysis': True,
    'do_NN_analysis': True
}
config = dict2namespace(params)

if config.do_sample_analysis:
    os.chdir(config.battery_path)

    if os.path.exists(config.results_df_path):
        results_df = pd.read_pickle(config.results_df_path)
    else:
        results_df = pd.DataFrame(columns=["run_num"])

    dirs = os.listdir()

    shuffle(dirs)

    for run_dir in tqdm.tqdm(dirs):  # loop over run directories in the battery
        os.chdir(config.battery_path)

        go_forward = False
        try:
            int(run_dir)
            go_forward = True
        except:
            pass

        if go_forward:

            os.chdir(run_dir)

            params_dict = np.load('run_config.npy', allow_pickle=True).item()

            new_row = {
                "temperature_series": [np.zeros(1)],
                "pressure_series": [np.zeros(1)],
                "E_pair": [np.zeros(1)],
                "E_mol": [np.zeros(1)],
                "E_tot": [np.zeros(1)],
                "intermolecular_rdfs": [np.zeros(1)],
                "rdf_drift": [np.zeros(1)],
                "rdf_times": [np.zeros(1)],
                "global_Ip_alignment": [np.zeros(1)],
                "local_Ip_alignment": [np.zeros(1)],
            }
            new_row.update(params_dict)

            # check if the run crashed
            # files = os.listdir()
            # for file in files:
            #     if 'slurm-' in file:
            #         slurm_filename = file
            #         break
            # reader = open(slurm_filename, 'r')
            # text = reader.read()
            # reader.close()

            if os.path.exists('new_traj.dump') and (int(run_dir) not in list(results_df['run_num'])):
                # do the analysis
                u = mda.Universe("system.data", "traj.dcd", format="LAMMPS")

                # find the reference system for comparison
                current_size = new_row['cluster_size']
                current_structure = new_row['structure_identifier']
                current_temperature = new_row['temperature']
                new_row['cluster_size'] = [new_row['cluster_size']]  # necessary for this to be a list inside a list

                print(run_dir)

                '''thermodynamic data'''
                thermo_results_dict = process_thermo_data()
                thermo_figs = plot_thermodynamic_data(thermo_results_dict)
                new_row["temperature_series"] = [thermo_results_dict['temp']]
                new_row["pressure_series"] = [thermo_results_dict['Press']]
                new_row["E_pair"] = [thermo_results_dict["E_pair"]]
                new_row["E_mol"] = [thermo_results_dict["E_mol"]]
                new_row["E_tot"] = [thermo_results_dict["E_tot"]]
                new_row["ns_per_day"] = [int(thermo_results_dict["ns_per_day"])]

                if config.do_alignment_analysis:
                    '''alignment analysis'''
                    Ip_trajectory, Ip_overlap_trajectory, alignment_frames = cluster_molecule_alignment(u, print_steps=25)  # len(u.trajectory))
                    global_molwise_alignment, local_molwise_alignment = compute_Ip_molwise_alignment(u, Ip_overlap_trajectory, alignment_frames)

                    if config.write_trajectory:
                        # convert from atomwise to molwise trajectory
                        atomwise_alignment_trajectory = np.concatenate([local_molwise_alignment[:, rind, None].repeat(u.residues[rind].atoms.n_atoms, 1) for rind in range(u.residues.n_residues)], axis=-1)
                        # local_molwise_alignment.repeat(15, 1)

                        rewrite_trajectory(u, run_dir, extra_atomwise_values=atomwise_alignment_trajectory)  # todo update repeat pattern for mixed benzamides

                    new_row["global_Ip_alignment"] = [global_molwise_alignment.mean((1, 2))]
                    new_row["local_Ip_alignment"] = [local_molwise_alignment.mean(1)]

                if config.do_rdf_analysis:
                    '''rdf analysis'''
                    atomwise_rdfs, bins, rdf_times = trajectory_rdf_analysis(
                        u, nbins=100, rrange=[0, 5], print_steps=25, core_cutoff=7.5)
                    new_row["intermolecular_rdfs"] = [atomwise_rdfs]
                    new_row["rdf_times"] = [rdf_times]
                    new_row["rdf_bins"] = [bins]

                if config.do_NN_analysis:
                    NN_output = process_dump('new_traj.dump')

                    '''compute following properties
                    time and mol wise averages
                    mol-averaged trajectories
                    mol-wise variability
                    '''
                    all_steps = pd.concat([frame for frame in NN_output.values()])

                    categories = [column for column in all_steps.columns if 'NNout' in column]
                    avg_probs = [np.mean(all_steps[column]) for column in all_steps.columns if 'NNout' in column]
                    prob_trajectories = np.asarray([[np.mean(df[column]) for column in df.columns if 'NNout' in column] for df in NN_output.values()])

                    n_mols = len(np.unique(list(NN_output.values())[0]['mol']))
                    horiz_combo = pd.concat([frame for frame in NN_output.values()], axis=1)

                    atomwise_variance = np.asarray([np.var(horiz_combo[column], axis=1) for column in categories])
                    #
                    # '''make figs'''
                    # fig = make_subplots(rows=1, cols=3, subplot_titles=['Trajectory Mean', 'Variability', 'Trajectory'])
                    # fig.add_trace(go.Bar(x=categories, y=avg_probs, showlegend=False),
                    #               row=1, col=1)
                    # fig.add_trace(go.Bar(x=categories, y=atomwise_variance.mean(1), showlegend=False),
                    #               row=1, col=2)
                    # for i in range(len(categories)):
                    #     fig.add_scattergl(y=prob_trajectories[:, i], name=categories[i]),
                    #                   row=1, col=3)
                    # fig.show()
                    #
                    # fig.write_image('NNout.png')

                    new_row['NN_classes'] = [categories]
                    new_row['NN_trajectories'] = [prob_trajectories]
                    new_row['NN_means'] = [avg_probs]
                    new_row['NN_variance'] = [atomwise_variance.mean(1)]

                if os.path.exists(config.results_df_path):
                    results_df = pd.read_pickle(config.results_df_path)

                results_df = pd.concat([results_df, pd.DataFrame.from_dict(new_row)], axis=0)
                results_df.to_pickle('../' + config.results_df_path)

aa = 1

if config.do_reference_analysis:
    os.chdir(config.reference_path)

    if os.path.exists(config.reference_df_path):
        results_df = pd.read_pickle(config.reference_df_path)
    else:
        results_df = pd.DataFrame(columns=["temperature"])  # todo need unique numbering system for runs

    dirs = os.listdir()

    for run_dir in dirs:  # loop over run directories in the battery
        if 'nvt' in run_dir:
            os.chdir(run_dir)
            if '13' in run_dir:
                structure_identifier = 'NICOAM13'
            elif '17' in run_dir:
                structure_identifier = "NICOAM17"
            elif '16' in run_dir:
                structure_identifier = "NICOAM16"

            files = os.listdir()
            for file in files:
                if 'dcd' in file and 'nvt' in file:
                    temperature = int(file.split('_')[0][1:])

                    new_row = {
                        "structure_identifier": [structure_identifier],
                        "temperature": [temperature],
                        "global_Ip_alignment": [np.zeros(1)],
                        "local_Ip_alignment": [np.zeros(1)],
                    }

                    u = mda.Universe("system.data", file, format="LAMMPS")

                    print(file)

                    '''alignment analysis'''
                    if config.do_alignment_analysis:
                        '''alignment analysis'''
                        Ip_trajectory, Ip_overlap_trajectory, alignment_frames = cluster_molecule_alignment(u, print_steps=100)  # len(u.trajectory))
                        global_molwise_alignment, local_molwise_alignment = compute_Ip_molwise_alignment(u, Ip_overlap_trajectory, alignment_frames)

                        new_row["global_Ip_alignment"] = [global_molwise_alignment.mean((1, 2))]
                        new_row["local_Ip_alignment"] = [local_molwise_alignment.mean(1)]

                    if config.do_rdf_analysis:
                        '''rdf analysis'''
                        atomwise_rdfs, bins, rdf_times = trajectory_rdf_analysis(
                            u, nbins=100, rrange=[0, 5], print_steps=25, core_cutoff=10)

                        new_row["intermolecular_rdfs"] = [atomwise_rdfs]
                        new_row["rdf_times"] = [rdf_times]

                    results_df = pd.concat([results_df, pd.DataFrame.from_dict(new_row)])

            os.chdir('../')
    results_df.to_pickle(config.reference_df_path)

aa = 1

import MDAnalysis as mda
import os
import wandb
from cluster_figs import \
    (plot_rdf_series, plot_intermolecular_rdf_series,
     plot_cluster_stability, plot_cluster_centroids_drift,
     process_thermo_data, plot_atomwise_rdf_drift, plot_alignment_fingerprint,
     plot_thermodynamic_data, trajectory_rdf_analysis,
     plot_atomwise_rdf_ref_dist)
from utils import (dict2namespace, names_dict, ff_names_dict)
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

params = {
    'reference_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\bulk_reference2/',
    'battery_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\battery_9/',
    'machine': 'local',  # or 'cluster'  ### doesn't do anything
    'show_figs': False,
    'write_trajectory': False,
    'make_run_wise_figs': True,
}
config = dict2namespace(params)

wandb.init(config=params, project="nicotinamide_clusters",
           entity="mkilgour", tags=["bulk_reference_test"],
           settings=wandb.Settings(code_dir="."))

wandb.run.name = config.battery_path
wandb.run.save()

os.chdir(config.battery_path)

if os.path.exists('results_df'):
    results_df = pd.read_pickle('results_df')
else:
    results_df = pd.DataFrame(columns=["run_num",
                                       "temperature",
                                       "pressure",
                                       "E_pair",
                                       "E_mol",
                                       "E_tot",
                                       "full_rdf",
                                       "full_intermolecular_rdf",
                                       "atomwise_intermolecular_rdfs",
                                       ])
'''reference indexes'''
temperatures = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
crystal_structures = ["NICOAM13", "NICOAM17"]

ref_temp_list = []
ref_crystal_list = []
for j in range(len(temperatures)):
    for k in range(len(crystal_structures)):
        ref_temp_list.append(temperatures[j])
        ref_crystal_list.append(crystal_structures[k])

'''run indexes'''
cluster_sizes = [[2, 2, 2],
                 [3, 3, 3],
                 [4, 4, 4],
                 [2, 2, 4],
                 [2, 4, 2],
                 [4, 2, 2]]
temperatures = [150, 300, 450]
crystal_structures = ["NICOAM13", "NICOAM17"]

n_runs = len(cluster_sizes) * len(temperatures) * len(crystal_structures)
run_nums = list(np.arange(1, n_runs + 1))

ind = 0
size_list = []
temp_list = []
crystal_list = []
for i in range(len(cluster_sizes)):
    for j in range(len(temperatures)):
        for k in range(len(crystal_structures)):
            size_list.append(cluster_sizes[i])
            temp_list.append(temperatures[j])
            crystal_list.append(crystal_structures[k])

dirs = os.listdir()
for run_dir in dirs:  # loop over run directories in the battery
    os.chdir(config.battery_path)

    if (run_dir != 'common') and \
            (run_dir not in results_df["run_num"].values) and \
            ('results_df' not in run_dir) and \
            ('png' not in run_dir):
        os.chdir(run_dir)
        # do the analysis
        u = mda.Universe("system.data", "traj.dcd", format="LAMMPS")
        # find the reference system for comparison
        current_structure = crystal_list[int(run_dir)-1]
        current_temperature = temp_list[int(run_dir)-1]

        # index for relevant reference system
        ref_index = np.argwhere((np.asarray(ref_temp_list) == current_temperature) * (np.asarray(ref_crystal_list) == current_structure))[0][0]
        ref_path = params['reference_path'] + str(ref_index + 1)
        ref_u = mda.Universe(ref_path + "/system.data", ref_path + "/traj.dcd", format="LAMMPS")

        print(run_dir)

        if config.write_trajectory:
            if not os.path.exists(f"{run_dir}_traj.xyz"):
                atom_types = u.atoms.types
                atom_names = np.asarray([names_dict[atype] for atype in atom_types])
                u.add_TopologyAttr('name', atom_names)
                cluster = u.select_atoms("all")
                with mda.Writer(f"{run_dir}_traj.xyz", cluster.n_atoms) as W:
                    for ts in u.trajectory:
                        W.write(cluster)

        if config.make_run_wise_figs:
            '''thermodynamic data'''
            thermo_results_dict = process_thermo_data()
            fig = plot_thermodynamic_data(thermo_results_dict)
            if config.show_figs:
                fig.show(renderer="browser")
            wandb.log({'Thermo Data': fig})
            fig.write_image('Thermo Data.png')

            '''rdf analysis'''
            full_rdf, intermolecular_rdf, atomwise_rdfs, bins = trajectory_rdf_analysis(u)
            '''reference rdf analysis'''
            ref_full_rdf, ref_intermolecular_rdf, ref_atomwise_rdfs, bins = trajectory_rdf_analysis(ref_u)

            '''intermolecular atomwise rdf distances'''
            fig, rdf_times, rdf_drift = plot_atomwise_rdf_ref_dist(u, atomwise_rdfs, ref_atomwise_rdfs, bins)
            if config.show_figs:
                fig.show(renderer="browser")
            wandb.log({'Intermolecular Atomwise RDF Drift': fig})
            fig.write_image('cluster_rdf_shift.png')

            '''save results'''
            new_row = {"run_num": run_dir,
                       "temperature": [thermo_results_dict['temp']],
                       "pressure": [thermo_results_dict['Press']],
                       "E_pair": [thermo_results_dict["E_pair"]],
                       "E_mol": [thermo_results_dict["E_mol"]],
                       "E_tot": [thermo_results_dict["E_tot"]],
                       "full_rdf": [full_rdf],
                       "full_intermolecular_rdf": [intermolecular_rdf],
                       "atomwise_intermolecular_rdfs": [atomwise_rdfs],
                       "rdf_drift": [rdf_drift],
                       "rdf_times": [rdf_times],
                       }
            results_df = pd.concat([results_df, pd.DataFrame.from_dict(new_row)])
            results_df.to_pickle('../results_df')

aa = 0


n_runs = 36
cluster_sizes = [[2, 2, 2],
                 [3, 3, 3],
                 [4, 4, 4],
                 [2, 2, 4],
                 [2, 4, 2],
                 [4, 2, 2]]
temperatures = [150, 300, 450]
crystal_structures = ["NICOAM13", "NICOAM17"]

n_runs = len(cluster_sizes) * len(temperatures) * len(crystal_structures)
run_nums = list(np.arange(1, n_runs + 1))

ind = 0
size_list = []
temp_list = []
crystal_list = []
for i in range(len(cluster_sizes)):
    for j in range(len(temperatures)):
        for k in range(len(crystal_structures)):
            size_list.append(cluster_sizes[i])
            temp_list.append(temperatures[j])
            crystal_list.append(crystal_structures[k])

results_df['run size'] = [size_list[int(val) - 1] for val in results_df['run_num'].values]
results_df['run crystal'] = [crystal_list[int(val) - 1] for val in results_df['run_num'].values]
results_df['run temperature'] = [temp_list[int(val) - 1] for val in results_df['run_num'].values]

'''RDF shift'''
shift_heatmap = np.zeros((len(crystal_structures), len(cluster_sizes), len(temperatures)))
for iT, temp in enumerate(temperatures):
    for iS, cs in enumerate(cluster_sizes):
        for iC, cr in enumerate(crystal_structures):
            for ii, row in results_df.iterrows():
                if row['run size'] == cs:
                    if row['run crystal'] == cr:
                        if row['run temperature'] == temp:
                            shift_heatmap[iC, iS, iT] = (row['rdf_drift'].mean())

fig = make_subplots(rows=1, cols=2, subplot_titles=crystal_structures)
fig.add_trace(go.Heatmap(z=np.log10(shift_heatmap[0])), row=1, col=1)
fig.add_trace(go.Heatmap(z=np.log10(shift_heatmap[1])), row=1, col=2)
fig.update_xaxes(title_text='Temperature', row=1, col=1)
fig.update_yaxes(title_text='Cluster Size', row=1, col=1)
fig.update_xaxes(title_text='Temperature', row=1, col=2)
fig.update_yaxes(title_text='Cluster Size', row=1, col=2)
fig.update_layout(xaxis=dict(
    tickmode='array',
    tickvals=[0, 1, 2],
    ticktext=temperatures
))
fig.update_layout(yaxis=dict(
    tickmode='array',
    tickvals=[0, 1, 2, 3, 4, 5],
    ticktext=cluster_sizes
))
fig.update_layout(xaxis2=dict(
    tickmode='array',
    tickvals=[0, 1, 2],
    ticktext=temperatures
))
fig.update_layout(yaxis2=dict(
    tickmode='array',
    tickvals=[0, 1, 2, 3, 4, 5],
    ticktext=cluster_sizes
))
fig.update_layout(title="log10 RDF Shift")
fig.show(renderer="browser")
fig.write_image("rdf_shift_heatmap.png")

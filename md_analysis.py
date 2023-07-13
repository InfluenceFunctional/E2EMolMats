import MDAnalysis as mda
import os
import wandb
from cluster_figs import \
    (plot_rdf_series, plot_intermolecular_rdf_series,
     plot_cluster_stability, plot_cluster_centroids_drift,
     process_thermo_data, plot_atomwise_rdf_drift, plot_alignment_fingerprint)
from utils import (dict2namespace, names_dict, ff_names_dict)
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd

params = {
    'battery_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\test_battery/',
    'machine': 'local',  # or 'cluster'  ### doesn't do anything
    'show_figs': False,
    'write_trajectory': False,
    'make_run_wise_figs': True,
}
config = dict2namespace(params)

wandb.init(config=params, project="nicotinamide_clusters",
           entity="mkilgour", tags=["reporting_test2"],
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
                                       "cluster_stability",
                                       "cluster_intermolecular_drift",
                                       "stability times",
                                       "distmat times",
                                       ])

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
            fig = make_subplots(rows=2, cols=3)
            ind = 0
            for i, key in enumerate(thermo_results_dict.keys()):
                if key != 'time step':
                    ind += 1
                    row = ind // 3 + 1
                    col = ind % 3 + 1
                    fig.add_trace(
                        go.Scattergl(x=thermo_results_dict['time step'],
                                     y=thermo_results_dict[key], name=key),
                        row=row, col=col
                    )
            if config.show_figs:
                fig.show()
            wandb.log({'Thermo Data': fig})
            fig.write_image('Thermo Data.png')

            '''full rdf'''
            atom_types = u.atoms.types
            atom_names = np.asarray([ff_names_dict[atype] for atype in atom_types])
            u.add_TopologyAttr('name', atom_names)
            fig, bins, full_rdf = plot_rdf_series(u)
            if config.show_figs:
                fig.show()
            wandb.log({'RDF Series': fig})
            fig.write_image('RDF_series.png')

            '''intermolecular df'''
            fig, bins, intermolecular_rdf = plot_intermolecular_rdf_series(u)
            if config.show_figs:
                fig.show()
            wandb.log({'Intermolecular RDF Series': fig})
            fig.write_image('intermolecular_RDF_series.png')

            '''atom-type-wise intermolecular rdf'''
            atomwise_rdfs = []
            for atom_type in ff_names_dict.values():
                if 'h' not in atom_type:
                    fig, bins, atom_type_rdf = plot_intermolecular_rdf_series(u, atom_type, atom_type, n_frames_avg=1)
                    if config.show_figs:
                        fig.show()
                    wandb.log({atom_type + ' Intermolecular RDF Series': fig})
                    fig.write_image(atom_type + ' intermolecular_RDF_series.png')
                    atomwise_rdfs.append(atom_type_rdf)

            '''cluster stability'''
            fig, stability_times, majority_proportion_fractions = plot_cluster_stability(u)
            if config.show_figs:
                fig.show()
            wandb.log({'Cluster Stability': fig})
            fig.write_image('cluster_stability.png')

            '''intramolecular centroids fingerprint'''
            fig, distmat_times, distmat_drift = plot_cluster_centroids_drift(u)
            if config.show_figs:
                fig.show()
            wandb.log({'Cluster Centroids Drift': fig})
            fig.write_image('cluster_centroids_drift.png')

            '''intermolecular atomwise rdf distances'''
            fig, rdf_times, rdf_drift = plot_atomwise_rdf_drift(u, atomwise_rdfs, bins)
            if config.show_figs:
                fig.show()
            wandb.log({'Intermolecular Atomwise RDF Drift': fig})
            fig.write_image('cluster_rdf_shift.png')

            '''molecule alignment fingerprint'''
            fig, Ip_overlaps_drift, Ip_overlaps_times = plot_alignment_fingerprint(u)
            if config.show_figs:
                fig.show()
            wandb.log({'Molecule Alignment Drift': fig})
            fig.write_image('cluster_alignment_shift.png')

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
                       "cluster_stability": [majority_proportion_fractions],
                       "cluster_intermolecular_drift": [distmat_drift],
                       "stability times": [stability_times],
                       "distmat times": [distmat_times],
                       "rdf_drift": [rdf_drift],
                       "rdf times": [rdf_times],
                       "Ip_overlaps_drift": [Ip_overlaps_drift],
                       "Ip overlaps times": [Ip_overlaps_times],
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

'''stability fig'''
stability_heatmap = np.zeros((len(crystal_structures), len(cluster_sizes), len(temperatures)))
for iT, temp in enumerate(temperatures):
    for iS, cs in enumerate(cluster_sizes):
        for iC, cr in enumerate(crystal_structures):
            for ii, row in results_df.iterrows():
                if row['run size'] == cs:
                    if row['run crystal'] == cr:
                        if row['run temperature'] == temp:
                            stability_heatmap[iC, iS, iT] = row['cluster_stability'][-3:].mean()

fig = make_subplots(rows=1, cols=2, subplot_titles=crystal_structures)
fig.add_trace(go.Heatmap(z=stability_heatmap[0]), row=1, col=1)
fig.add_trace(go.Heatmap(z=stability_heatmap[1]), row=1, col=2)
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
fig.update_layout(title="Cluster Stability")
fig.show()
fig.write_image("overall_cluster_stability.png")

'''distmat shift'''
shift_heatmap = np.zeros((len(crystal_structures), len(cluster_sizes), len(temperatures)))
for iT, temp in enumerate(temperatures):
    for iS, cs in enumerate(cluster_sizes):
        for iC, cr in enumerate(crystal_structures):
            for ii, row in results_df.iterrows():
                if row['run size'] == cs:
                    if row['run crystal'] == cr:
                        if row['run temperature'] == temp:
                            shift_heatmap[iC, iS, iT] = np.log10(row['cluster_intermolecular_drift'].mean())

fig = make_subplots(rows=1, cols=2, subplot_titles=crystal_structures)
fig.add_trace(go.Heatmap(z=shift_heatmap[0]), row=1, col=1)
fig.add_trace(go.Heatmap(z=shift_heatmap[1]), row=1, col=2)
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
fig.update_layout(title="log Intermolecular Distances Shift")
fig.show()
fig.write_image("overall_distances_shift.png")


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
fig.add_trace(go.Heatmap(z=shift_heatmap[0]), row=1, col=1)
fig.add_trace(go.Heatmap(z=shift_heatmap[1]), row=1, col=2)
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
fig.update_layout(title="RDF Shift")
fig.show()
fig.write_image("rdf_shift_heatmap.png")

'''Ip_overlaps'''
shift_heatmap = np.zeros((len(crystal_structures), len(cluster_sizes), len(temperatures)))
for iT, temp in enumerate(temperatures):
    for iS, cs in enumerate(cluster_sizes):
        for iC, cr in enumerate(crystal_structures):
            for ii, row in results_df.iterrows():
                if row['run size'] == cs:
                    if row['run crystal'] == cr:
                        if row['run temperature'] == temp:
                            shift_heatmap[iC, iS, iT] = (row['Ip_overlaps_drift'].mean())

fig = make_subplots(rows=1, cols=2, subplot_titles=crystal_structures)
fig.add_trace(go.Heatmap(z=shift_heatmap[0]), row=1, col=1)
fig.add_trace(go.Heatmap(z=shift_heatmap[1]), row=1, col=2)
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
fig.update_layout(title="Ip Overlap Shift")
fig.show()
fig.write_image("Ip_overlaps_shift_heatmap.png")

import MDAnalysis as mda
import os
import wandb

from analysis.residue_test import residue_cleanup
from reporting.cluster_figs import \
    (plot_rdf_series, plot_intermolecular_rdf_series,
     plot_cluster_stability, plot_cluster_centroids_drift,
     process_thermo_data, plot_atomwise_rdf_drift, plot_alignment_fingerprint,
     plot_thermodynamic_data, trajectory_rdf_analysis,
     plot_atomwise_rdf_ref_dist, cluster_molecule_alignment)
from utils import (dict2namespace, names_dict, ff_names_dict, cell_vol, rewrite_trajectory, tile_universe, compute_Ip_alignment)
import numpy as np
from plotly.subplots import make_subplots
from scipy.spatial.distance import cdist, pdist
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import plotly.io as pio

pio.renderers.default = 'browser'

params = {
    'reference_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\bulk_reference4/',
    'battery_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\benzamide_test2/',
    'machine': 'local',  # or 'cluster'  ### doesn't do anything
    'show_figs': False,
    'write_trajectory': False,
    'make_run_wise_figs': True,
}
config = dict2namespace(params)

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

'''run indices'''

cluster_sizes = [[4, 4, 4]]
temperatures = [250, 300, 350]
crystal_structures = ["NICOAM13", "NICOAM17"]
defect_rates = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5]

n_runs = len(cluster_sizes) * len(temperatures) * len(crystal_structures) * len(defect_rates)
run_nums = list(np.arange(1, n_runs + 1))

ind = 0
size_list = []
temp_list = []
crystal_list = []
defect_list = []
for i in range(len(cluster_sizes)):
    for j in range(len(temperatures)):
        for k in range(len(crystal_structures)):
            for l in range(len(defect_rates)):
                size_list.append(cluster_sizes[i])
                temp_list.append(temperatures[j])
                crystal_list.append(crystal_structures[k])
                defect_list.append(defect_rates[l])

params_dict = {
    'defect_rate': defect_list,
    'size': size_list,
    'temperature': temp_list,
    'crystal': crystal_list,
}

params_df = pd.DataFrame.from_dict(params_dict)

dirs = os.listdir()

# wandb.init(config=params, project="nicotinamide_clusters",
#            entity="mkilgour", tags=["bulk_reference_test"],
#            settings=wandb.Settings(code_dir="."))
#
# wandb.run.name = config.battery_path
# wandb.run.save()

for run_dir in dirs:  # loop over run directories in the battery
    os.chdir(config.battery_path)

    if (run_dir != 'common') and \
            (run_dir not in results_df["run_num"].values) and \
            ('results_df' not in run_dir) and \
            ('png' not in run_dir) and \
            ('wandb' not in run_dir):
        os.chdir(run_dir)

        # check if the run crashed
        files = os.listdir()
        for file in files:
            if 'slurm-' in file:
                slurm_filename = file
                break
        reader = open(slurm_filename, 'r')
        text = reader.read()
        reader.close()

        if 'oom' in text:  # run crashed
            '''save results'''
            new_row = {"run_num": run_dir,
                       "temperature": [np.zeros(1)],
                       "pressure": [np.zeros(1)],
                       "E_pair": [np.zeros(1)],
                       "E_mol": [np.zeros(1)],
                       "E_tot": [np.zeros(1)],
                       "atomwise_intermolecular_rdfs": [np.zeros(1)],
                       "rdf_drift": [np.zeros(1)],
                       "rdf_times": [np.zeros(1)],
                       }
        else:
            # do the analysis

            # sort out residues manually
            if 'benzamide_test2' in params['battery_path']:
                residue_cleanup()
                u = mda.Universe("fixed_system.data", "traj.dcd", format="LAMMPS")
            else:
                #automated residues alignment
                u = mda.Universe("system.data", "traj.dcd", format="LAMMPS")


            # find the reference system for comparison
            current_size = size_list[int(run_dir) - 1]
            current_structure = crystal_list[int(run_dir) - 1]
            current_temperature = temp_list[int(run_dir) - 1]

            # index for relevant reference system
            ref_index = np.argwhere((np.asarray(ref_temp_list) == current_temperature) * (np.asarray(ref_crystal_list) == current_structure))[0][0]
            ref_path = params['reference_path'] + str(ref_index + 1)
            ref_u = mda.Universe(ref_path + "/system.data", ref_path + "/traj.dcd", format="LAMMPS")

            # get reference and sample density
            coords = u.atoms.positions
            coords -= coords.mean(0)
            dists = cdist(np.asarray((0, 0, 0))[None, :], coords)
            subbox_size = min(10, u.dimensions[:3].min() / 4)

            density = np.sum(dists < subbox_size) / ((4 / 3) * np.pi * subbox_size ** 3)
            ref_density = len(ref_u.atoms) / cell_vol(ref_u.dimensions[:3], ref_u.dimensions[3:], units='degrees')
            density_difference = np.abs((ref_density - density)) / ref_density

            print(run_dir)
            if config.write_trajectory:
                rewrite_trajectory(u, run_dir)

            '''thermodynamic data'''
            thermo_results_dict = process_thermo_data()
            thermo_fig = plot_thermodynamic_data(thermo_results_dict)

            '''alignment analysis'''
            Ip_trajectory, Ip_overlap_trajectory = cluster_molecule_alignment(u)
            ref_Ip_trajectory, ref_Ip_overlap_trajectory = cluster_molecule_alignment(ref_u)

            Ip_alignment_trajectory = compute_Ip_alignment(u, Ip_overlap_trajectory)
            ref_Ip_alignment_trajectory =compute_Ip_alignment(ref_u, ref_Ip_overlap_trajectory)

            '''reference rdf analysis'''
            ref_atomwise_rdfs, bins, rdf_times = trajectory_rdf_analysis(
                ref_u, nbins=100, rrange=[0, 8], core_cutoff=subbox_size, tiling=current_size)

            '''rdf analysis'''
            atomwise_rdfs, bins, rdf_times = trajectory_rdf_analysis(
                u, nbins=100, rrange=[0, 8], core_cutoff=subbox_size)

            '''intermolecular atomwise rdf distances'''
            rdf_drift = plot_atomwise_rdf_ref_dist(u, atomwise_rdfs, ref_atomwise_rdfs, bins)

            '''save results'''
            new_row = {"run_num": run_dir,
                       "temperature": [thermo_results_dict['temp']],
                       "pressure": [thermo_results_dict['Press']],
                       "E_pair": [thermo_results_dict["E_pair"]],
                       "E_mol": [thermo_results_dict["E_mol"]],
                       "E_tot": [thermo_results_dict["E_tot"]],
                       "ns/day": [thermo_results_dict["ns/day"]],
                       "atomwise_intermolecular_rdfs": [atomwise_rdfs],
                       "rdf_drift": [rdf_drift],
                       "rdf_times": [rdf_times],
                       #"Ip_trajectory": [Ip_trajectory],
                       "Ip_alignment": [Ip_alignment_trajectory],
                       "ref_Ip_alignment": [ref_Ip_alignment_trajectory],
                       #"ref_Ip_trajectory": [ref_Ip_trajectory],
                       }
        results_df = pd.concat([results_df, pd.DataFrame.from_dict(new_row)])
        results_df.to_pickle('../results_df')
        #
        # if config.show_figs:
        #     thermo_fig.show(renderer="browser")
        # wandb.log({'Thermo Data': thermo_fig})
        # thermo_fig.write_image('Thermo Data.png')

aa = 0

results_df['run size'] = [size_list[int(val) - 1] for val in results_df['run_num'].values]
results_df['run crystal'] = [crystal_list[int(val) - 1] for val in results_df['run_num'].values]
results_df['run temperature'] = [temp_list[int(val) - 1] for val in results_df['run_num'].values]
results_df['run vacancy rate'] = [defect_list[int(val) - 1] for val in results_df['run_num'].values]

if "Ip_alignment" in results_df.columns:
    '''
    whether inertial axes are aligned over time - in Z'=4 structure perfect alignment will be 25%
    '''
    shift_heatmap = np.zeros((len(crystal_structures), len(defect_rates), len(temperatures)))
    for iT, temp in enumerate(temperatures):
        for iC, cr in enumerate(crystal_structures):
            for iG, gr in enumerate(defect_rates):
                for ii, row in results_df.iterrows():
                    if row['run crystal'] == cr:
                        if row['run temperature'] == temp:
                            if row['run vacancy rate'] == gr:
                                try:
                                    shift_heatmap[iC, iG, iT] = row['Ip_alignment'].mean()# / row['ref_Ip_alignment'].mean()
                                except:
                                    shift_heatmap[iC, iG, iT] = 0

    fig = make_subplots(rows=1, cols=2, subplot_titles=crystal_structures)
    fig.add_trace(go.Heatmap(z=(shift_heatmap[0])), row=1, col=1)
    fig.add_trace(go.Heatmap(z=(shift_heatmap[1])), row=1, col=2)
    fig.update_xaxes(title_text='Temperature', row=1, col=1)
    fig.update_yaxes(title_text='defect_rate', row=1, col=1)
    fig.update_xaxes(title_text='Temperature', row=1, col=2)
    fig.update_yaxes(title_text='defect_rate', row=1, col=2)
    fig.update_layout(xaxis=dict(
        tickmode='array',
        tickvals=np.arange(len(temperatures)),
        ticktext=temperatures
    ))
    fig.update_layout(yaxis=dict(
        tickmode='array',
        tickvals=np.arange(len(defect_rates)),
        ticktext=defect_rates
    ))
    fig.update_layout(xaxis2=dict(
        tickmode='array',
        tickvals=np.arange(len(temperatures)),
        ticktext=temperatures
    ))
    fig.update_layout(yaxis2=dict(
        tickmode='array',
        tickvals=np.arange(len(defect_rates)),
        ticktext=defect_rates
    ))
    fig.update_layout(title="Cluster Alignment")
    fig.show(renderer="browser")
    fig.write_image("alignment_shift_heatmap.png")

if "rdf_drift" in results_df.columns:
    '''RDF shift'''
    shift_heatmap = np.zeros((len(crystal_structures), len(defect_rates), len(temperatures)))
    for iT, temp in enumerate(temperatures):
        for iC, cr in enumerate(crystal_structures):
            for iG, gr in enumerate(defect_rates):
                for ii, row in results_df.iterrows():
                    if row['run crystal'] == cr:
                        if row['run temperature'] == temp:
                            if row['run vacancy rate'] == gr:
                                try:
                                    shift_heatmap[iC, iG, iT] = (row['rdf_drift'].mean())
                                except:
                                    shift_heatmap[iC, iG, iT] = 0

    fig = make_subplots(rows=1, cols=2, subplot_titles=crystal_structures)
    fig.add_trace(go.Heatmap(z=(shift_heatmap[0])), row=1, col=1)
    fig.add_trace(go.Heatmap(z=(shift_heatmap[1])), row=1, col=2)
    fig.update_xaxes(title_text='Temperature', row=1, col=1)
    fig.update_yaxes(title_text='defect_rate', row=1, col=1)
    fig.update_xaxes(title_text='Temperature', row=1, col=2)
    fig.update_yaxes(title_text='defect_rate', row=1, col=2)
    fig.update_layout(xaxis=dict(
        tickmode='array',
        tickvals=np.arange(len(temperatures)),
        ticktext=temperatures
    ))
    fig.update_layout(yaxis=dict(
        tickmode='array',
        tickvals=np.arange(len(defect_rates)),
        ticktext=defect_rates
    ))
    fig.update_layout(xaxis2=dict(
        tickmode='array',
        tickvals=np.arange(len(temperatures)),
        ticktext=temperatures
    ))
    fig.update_layout(yaxis2=dict(
        tickmode='array',
        tickvals=np.arange(len(defect_rates)),
        ticktext=defect_rates
    ))
    fig.update_layout(title="RDF Shift")
    fig.show(renderer="browser")
    fig.write_image("rdf_shift_heatmap.png")


if "ns/day" in results_df.columns:
    '''RDF shift'''
    shift_heatmap = np.zeros((len(crystal_structures), len(defect_rates), len(temperatures)))
    for iT, temp in enumerate(temperatures):
        for iC, cr in enumerate(crystal_structures):
            for iG, gr in enumerate(defect_rates):
                for ii, row in results_df.iterrows():
                    if row['run crystal'] == cr:
                        if row['run temperature'] == temp:
                            if row['run vacancy rate'] == gr:
                                try:
                                    shift_heatmap[iC, iG, iT] = (row['ns/day'])
                                except:
                                    shift_heatmap[iC, iG, iT] = 0

    fig = make_subplots(rows=1, cols=2, subplot_titles=crystal_structures)
    fig.add_trace(go.Heatmap(z=(shift_heatmap[0])), row=1, col=1)
    fig.add_trace(go.Heatmap(z=(shift_heatmap[1])), row=1, col=2)
    fig.update_xaxes(title_text='Temperature', row=1, col=1)
    fig.update_yaxes(title_text='defect_rate', row=1, col=1)
    fig.update_xaxes(title_text='Temperature', row=1, col=2)
    fig.update_yaxes(title_text='defect_rate', row=1, col=2)
    fig.update_layout(xaxis=dict(
        tickmode='array',
        tickvals=np.arange(len(temperatures)),
        ticktext=temperatures
    ))
    fig.update_layout(yaxis=dict(
        tickmode='array',
        tickvals=np.arange(len(defect_rates)),
        ticktext=defect_rates
    ))
    fig.update_layout(xaxis2=dict(
        tickmode='array',
        tickvals=np.arange(len(temperatures)),
        ticktext=temperatures
    ))
    fig.update_layout(yaxis2=dict(
        tickmode='array',
        tickvals=np.arange(len(defect_rates)),
        ticktext=defect_rates
    ))
    fig.update_layout(title="ns per day")
    fig.show(renderer="browser")
    fig.write_image("ns_day_heatmap.png")

aa = 1

# # sort out residues
# mol_type_ind = np.zeros(len(u.atoms),dtype=np.int64)
# atom_counter = 0
# mol_counter = 0
#
# for ind in range(1, len(u.atoms)):
#     atom_counter += 1
#     if (u.atoms[ind].type == '5' and u.atoms[ind-1].type == '1'):  # beginning of new molecule - assign indices to previous
#         if atom_counter == 16:  # benzamide
#             mol_type_ind[ind-16:ind] = mol_counter
#         elif atom_counter == 15:  # nicotinamide
#             mol_type_ind[ind-15:ind] = mol_counter
#
#         mol_counter += 1
#         atom_counter = 0
#
# # must also add final molecule
# last_index = np.max(np.argwhere(mol_type_ind == mol_counter - 1)[:, 0])
# last_mol_len = len(u.atoms) - last_index - 1
# if last_mol_len == 16:  # benzamide
#     mol_type_ind[-16:] = mol_counter
# elif last_mol_len == 15:  # nicotinamide
#     mol_type_ind[-15:] = mol_counter
#
# mol_groups = []  # collect residues
# for ind in range(max(mol_type_ind) + 1):
#     mol_groups.append(mda.AtomGroup(u.atoms[mol_type_ind == ind]))
#
# from MDAnalysis.core.universe import Merge
# u = Merge(*mol_groups)  # universe with correct molecule numbering

import MDAnalysis as mda
import os
import wandb

from reporting.cluster_figs import (
    plot_thermodynamic_data, trajectory_rdf_analysis,
    plot_atomwise_rdf_ref_dist, cluster_molecule_alignment, cluster_property_heatmap, process_thermo_data)
from utils import (dict2namespace, cell_vol, rewrite_trajectory, compute_Ip_alignment, compute_Ip_molwise_alignment)
import numpy as np
from plotly.subplots import make_subplots
from scipy.spatial.distance import cdist, pdist
import plotly.graph_objects as go
import pandas as pd
import plotly.io as pio
from scipy.ndimage import gaussian_filter1d
import tqdm

pio.renderers.default = 'browser'

params = {
    'reference_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\bulk_reference/',
    'battery_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\defect_clusters_1/',
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
    results_df = pd.DataFrame(columns=["run_num"])  # todo need unique numbering system for runs

'''reference indexes'''
temperatures = [10, 100, 200, 300, 350, 400, 450, 500, 550, 600, 700]
crystal_structures = ["NICOAM13", "NICOAM17"]

dirs = os.listdir()

for run_dir in tqdm.tqdm(dirs):  # loop over run directories in the battery
    os.chdir(config.battery_path)

    if (run_dir != 'common') and \
            (run_dir not in results_df["run_num"].values.astype(str)) and \
            ('results_df' not in run_dir) and \
            ('png' not in run_dir) and \
            ('wandb' not in run_dir):

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
        files = os.listdir()
        for file in files:
            if 'slurm-' in file:
                slurm_filename = file
                break
        reader = open(slurm_filename, 'r')
        text = reader.read()
        reader.close()

        if 'oom' not in text:
            # do the analysis
            u = mda.Universe("system.data", "traj.dcd", format="LAMMPS")

            # find the reference system for comparison
            current_size = new_row['cluster_size']
            current_structure = new_row['structure_identifier']
            current_temperature = new_row['temperature']

            print(run_dir)

            '''thermodynamic data'''
            thermo_results_dict = process_thermo_data()
            thermo_figs = plot_thermodynamic_data(thermo_results_dict)

            '''alignment analysis'''
            Ip_trajectory, Ip_overlap_trajectory, alignment_frames = cluster_molecule_alignment(u, print_steps=100) #len(u.trajectory))
            global_molwise_alignment, local_molwise_alignment = compute_Ip_molwise_alignment(u, Ip_overlap_trajectory, alignment_frames)

            if config.write_trajectory:
                # convert from atomwise to molwise trajectory
                atomwise_alignment_trajectory = np.concatenate([local_molwise_alignment[:,rind,None].repeat(u.residues[rind].atoms.n_atoms,1) for rind in range(u.residues.n_residues)], axis=-1)
                #local_molwise_alignment.repeat(15, 1)

                rewrite_trajectory(u, run_dir, extra_atomwise_values=atomwise_alignment_trajectory)# todo update repeat pattern for mixed benzamides

            #Ip_alignment_trajectory = compute_Ip_alignment(u, Ip_overlap_trajectory)
            #ref_Ip_alignment_trajectory = compute_Ip_alignment(ref_u, ref_Ip_overlap_trajectory)

            # '''reference rdf analysis'''
            # ref_atomwise_rdfs, bins, rdf_times = trajectory_rdf_analysis(
            #     ref_u, nbins=100, rrange=[0, 8], core_cutoff=subbox_size, tiling=current_size, print_steps=25)
            #
            # '''rdf analysis'''
            # atomwise_rdfs, bins, rdf_times = trajectory_rdf_analysis(
            #     u, nbins=100, rrange=[0, 8], core_cutoff=subbox_size, print_steps=25)
            #
            # '''intermolecular atomwise rdf distances'''
            # rdf_drift = plot_atomwise_rdf_ref_dist(u, atomwise_rdfs, ref_atomwise_rdfs, bins)

            '''save results'''
            new_row["temperature_series"] = [thermo_results_dict['temp']]
            new_row["pressure_series"] = [thermo_results_dict['Press']]
            new_row["E_pair"] = [thermo_results_dict["E_pair"]]
            new_row["E_mol"] = [thermo_results_dict["E_mol"]]
            new_row["E_tot"] = [thermo_results_dict["E_tot"]]
            new_row["ns_per_day"] = [thermo_results_dict["ns_per_day"]]
            #new_row["intermolecular_rdfs"] = [atomwise_rdfs]
            #new_row["rdf_drift"] = [rdf_drift]
            #new_row["rdf_times"] = [rdf_times]
            new_row["global_Ip_alignment"] = [global_molwise_alignment.mean((1,2))]
            new_row["local_Ip_alignment"] = [local_molwise_alignment.mean(1)]
            #new_row["normed_Ip_alignment"] = [Ip_alignment_trajectory / ref_Ip_alignment_trajectory.mean(0)]  # norm against average of ref timeseries
            new_row['cluster_size'] = [new_row['cluster_size']]

        results_df = pd.concat([results_df, pd.DataFrame.from_dict(new_row)])
        results_df.to_pickle('../results_df')

aa = 1
cluster_property_heatmap(results_df, 'global_Ip_alignment', 'temperature', 'defect_rate', take_mean=True)


#
# os.chdir(config.reference_path)
#
# if os.path.exists('reference_df'):
#     results_df = pd.read_pickle('reference_df')
# else:
#     results_df = pd.DataFrame(columns=["temperature"])  # todo need unique numbering system for runs
#
# dirs = os.listdir()
#
# for run_dir in dirs:  # loop over run directories in the battery
#     os.chdir(run_dir)
#     if '13' in run_dir:
#         structure_identifier = 'NICOAM13'
#     elif '17' in run_dir:
#         structure_identifier = "NICOAM17"
#
#     files = os.listdir()
#     for file in files:
#         if 'dcd' in file and 'nvt' in file:
#             temperature = int(file.split('_')[0][1:])
#
#             new_row = {
#                 "structure_identifier": [structure_identifier],
#                 "temperature": [temperature],
#                 "global_Ip_alignment": [np.zeros(1)],
#                 "local_Ip_alignment": [np.zeros(1)],
#             }
#
#             u = mda.Universe("system.data", file, format="LAMMPS")
#
#             print(file)
#
#             '''alignment analysis'''
#             Ip_trajectory, Ip_overlap_trajectory = cluster_molecule_alignment(u, print_steps=len(u.trajectory))
#             global_molwise_alignment, local_molwise_alignment = compute_Ip_molwise_alignment(u, Ip_overlap_trajectory)
#
#             new_row["global_Ip_alignment"] = [global_molwise_alignment.mean((1, 2))]
#             new_row["local_Ip_alignment"] = [local_molwise_alignment.mean(1)]
#
#             results_df = pd.concat([results_df, pd.DataFrame.from_dict(new_row)])
#             results_df.to_pickle('../reference_df')
#
#     os.chdir('../')
#
# aa = 1

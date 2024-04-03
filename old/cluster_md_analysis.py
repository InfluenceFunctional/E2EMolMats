import MDAnalysis as mda
import os

from e2emolmats.reporting import (
    plot_thermodynamic_data, plot_atomwise_rdf_ref_dist, cluster_molecule_alignment, cluster_property_heatmap)
from e2emolmats.reporting import process_thermo_data
from e2emolmats.reporting import trajectory_rdf_analysis
from e2emolmats.common.utils import (dict2namespace, rewrite_trajectory, compute_Ip_alignment, compute_Ip_molwise_alignment)
import numpy as np
from scipy.spatial.distance import cdist
import pandas as pd
import plotly.io as pio
from scipy.ndimage import gaussian_filter1d

pio.renderers.default = 'browser'

params = {
    'reference_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\bulk_reference/',
    'battery_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\dev8/',
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
temperatures = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
crystal_structures = ["NICOAM13", "NICOAM17"]

ref_temp_list = []
ref_crystal_list = []
for j in range(len(temperatures)):
    for k in range(len(crystal_structures)):
        ref_temp_list.append(temperatures[j])
        ref_crystal_list.append(crystal_structures[k])

dirs = os.listdir()

# wandb.init(config=params, project="E2EMolMats",
#            entity="mkilgour", tags=["bulk_reference_test"],
#            settings=wandb.Settings(code_dir="."))
#
# wandb.run.name = config.battery_path
# wandb.run.save()

for run_dir in dirs:  # loop over run directories in the battery
    os.chdir(config.battery_path)

    if (run_dir != 'md_data') and \
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
            "Ip_alignment": [np.zeros(1)],
            "ref_Ip_alignment": [np.zeros(1)],
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

        if not 'oom' in text:
            # do the analysis
            u = mda.Universe("system.data", "traj.dcd", format="LAMMPS")

            # find the reference system for comparison
            current_size = new_row['cluster_size']
            current_structure = new_row['structure_identifier']
            current_temperature = new_row['temperature']

            # index for relevant reference system
            # ref_index = np.argwhere((np.asarray(ref_temp_list) == current_temperature) * (np.asarray(ref_crystal_list) == current_structure))[0][0]
            # ref_path = params['reference_path'] + str(ref_index + 1)
            # ref_u = mda.Universe(ref_path + "/system.data", ref_path + "/traj.dcd", format="LAMMPS")

            # get reference and sample density
            coords = u.atoms.positions
            coords -= coords.mean(0)
            dists = cdist(np.asarray((0, 0, 0))[None, :], coords)
            subbox_size = min(10, np.ptp(u.atoms.positions) / 4)  # set the 'inside' molecules for which rdfs will be calculated

            # density = np.sum(dists < subbox_size) / ((4 / 3) * np.pi * subbox_size ** 3)
            # ref_density = len(ref_u.atoms) / cell_vol(ref_u.dimensions[:3], ref_u.dimensions[3:], units='degrees')
            # density_difference = np.abs((ref_density - density)) / ref_density

            print(run_dir)

            '''thermodynamic data'''
            thermo_results_dict = process_thermo_data()
            thermo_fig = plot_thermodynamic_data(thermo_results_dict)

            '''alignment analysis'''
            Ip_trajectory, Ip_overlap_trajectory = cluster_molecule_alignment(u, print_steps=len(u.trajectory))
            ref_Ip_trajectory, ref_Ip_overlap_trajectory = cluster_molecule_alignment(ref_u, print_steps=len(ref_u.trajectory))
            molwise_Ip_trajectory = compute_Ip_molwise_alignment(u, Ip_overlap_trajectory)

            smoothed_Ip_trajectory = gaussian_filter1d(molwise_Ip_trajectory.mean(-1), 5, axis=-1)
            extra_atomwise_values = smoothed_Ip_trajectory.repeat(15, 1)  # todo update repeat pattern for mixed benzamides

            if config.write_trajectory:
                rewrite_trajectory(u, run_dir, extra_atomwise_values=extra_atomwise_values)

            Ip_alignment_trajectory = compute_Ip_alignment(u, Ip_overlap_trajectory)
            ref_Ip_alignment_trajectory = compute_Ip_alignment(ref_u, ref_Ip_overlap_trajectory)

            '''reference rdf analysis'''
            ref_atomwise_rdfs, bins, rdf_times = trajectory_rdf_analysis(
                ref_u, nbins=100, rrange=[0, 8], core_cutoff=subbox_size, tiling=current_size, print_steps=25)

            '''rdf analysis'''
            atomwise_rdfs, bins, rdf_times = trajectory_rdf_analysis(
                u, nbins=100, rrange=[0, 8], core_cutoff=subbox_size, print_steps=25)

            '''intermolecular atomwise rdf distances'''
            rdf_drift = plot_atomwise_rdf_ref_dist(u, atomwise_rdfs, ref_atomwise_rdfs, bins)

            '''save results'''
            new_row["temperature_series"] = [thermo_results_dict['temp']]
            new_row["pressure_series"] = [thermo_results_dict['Press']]
            new_row["E_pair"] = [thermo_results_dict["E_pair"]]
            new_row["E_mol"] = [thermo_results_dict["E_mol"]]
            new_row["E_tot"] = [thermo_results_dict["E_tot"]]
            new_row["ns_per_day"] = [thermo_results_dict["ns_per_day"]]
            new_row["intermolecular_rdfs"] = [atomwise_rdfs]
            new_row["rdf_drift"] = [rdf_drift]
            new_row["rdf_times"] = [rdf_times]
            new_row["Ip_alignment"] = [Ip_alignment_trajectory]
            new_row["ref_Ip_alignment"] = [ref_Ip_alignment_trajectory]
            new_row["normed_Ip_alignment"] = [Ip_alignment_trajectory / ref_Ip_alignment_trajectory.mean(0)]  # norm against average of ref timeseries
            new_row['cluster_size'] = [new_row['cluster_size']]

        results_df = pd.concat([results_df, pd.DataFrame.from_dict(new_row)])
        results_df.to_pickle('../results_df')

aa = 0

results_df = results_df.reset_index()
cluster_property_heatmap(results_df, 'rdf_drift', 'temperature', 'defect_rate', take_mean=True)
cluster_property_heatmap(results_df, 'Ip_alignment', 'temperature', 'defect_rate', take_mean=True)
cluster_property_heatmap(results_df, 'ns_per_day', 'temperature', 'defect_rate', take_mean=False)

aa = 1

# # sort out residues

# # sort out residues manually
# if 'benzamide_test2' in params['battery_path']:
#     residue_cleanup()
#     u = mda.Universe("fixed_system.data", "traj.dcd", format="LAMMPS")
# else:
#     # automated residues alignment


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

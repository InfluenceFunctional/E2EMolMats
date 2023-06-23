import MDAnalysis as mda
import os
import wandb
from cluster_figs import \
    (plot_rdf_series, plot_intermolecular_rdf_series,
     plot_cluster_stability, plot_cluster_centroids_drift, process_thermo_data)
from utils import (dict2namespace, names_dict, ff_names_dict)
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go

params = {
    'battery_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\battery_8_1/',
    'machine': 'local',  # or 'cluster'  ### doesn't do anything
    'show_figs': False,
    'write_trajectory': True,
    'make_figs': True,
}
config = dict2namespace(params)

with wandb.init(config=params, project="nicotinamide_clusters",
                entity="mkilgour", tags=["reporting_test2"],
                settings=wandb.Settings(code_dir=".")):

    wandb.run.name = config.battery_path
    wandb.run.save()

    os.chdir(config.battery_path)
    dirs = os.listdir()
    for run_dir in dirs:  # loop over run directories in the battery
        os.chdir(config.battery_path)
        os.chdir(run_dir)
        # do the analysis
        u = mda.Universe("system.data", "traj.dcd", format="LAMMPS")

        if run_dir != 'common':
            if config.write_trajectory:
                atom_types = u.atoms.types

                atom_names = np.asarray([names_dict[atype] for atype in atom_types])
                u.add_TopologyAttr('name', atom_names)
                cluster = u.select_atoms("all")
                with mda.Writer(f"{run_dir}_traj.xyz", cluster.n_atoms) as W:
                    for ts in u.trajectory:
                        W.write(cluster)

            if config.make_figs:
                results_dict = process_thermo_data()
                fig = make_subplots(rows=2, cols=3)
                ind = 0
                for i, key in enumerate(results_dict.keys()):
                    if key != 'time step':
                        ind += 1
                        row = ind // 3 + 1
                        col = ind % 3 + 1
                        fig.add_trace(
                            go.Scattergl(x=results_dict['time step'],
                                         y=results_dict[key],name=key),
                            row=row, col=col
                        )
                if config.show_figs:
                    fig.show()
                wandb.log({'Thermo Data': fig})
                fig.write_image('Thermo Data.png')

                atom_types = u.atoms.types
                atom_names = np.asarray([ff_names_dict[atype] for atype in atom_types])
                u.add_TopologyAttr('name', atom_names)

                fig = plot_rdf_series(u)
                if config.show_figs:
                    fig.show()
                wandb.log({'RDF Series': fig})
                fig.write_image('RDF_series.png')

                fig = plot_intermolecular_rdf_series(u)
                if config.show_figs:
                    fig.show()
                wandb.log({'Intermolecular RDF Series': fig})
                fig.write_image('intermolecular_RDF_series.png')

                for atom_type in ff_names_dict.values():
                    fig = plot_intermolecular_rdf_series(u, atom_type, atom_type, n_frames_avg=10)
                    if config.show_figs:
                        fig.show()
                    wandb.log({atom_type + ' Intermolecular RDF Series': fig})
                    fig.write_image(atom_type + ' intermolecular_RDF_series.png')

                fig = plot_cluster_stability(u)
                if config.show_figs:
                    fig.show()
                wandb.log({'Cluster Stability': fig})
                fig.write_image('cluster_stability.png')

                fig = plot_cluster_centroids_drift(u)
                if config.show_figs:
                    fig.show()
                wandb.log({'Cluster Centroids Drift': fig})
                fig.write_image('cluster_centroids_drift.png')

                #  todo use mikes Ip axes and monitor dot product overlaps



    aa = 0

# cluster = u.select_atoms("all")
# with mda.Writer("traj.xyz", cluster.n_atoms) as W:
#     for ts in u.trajectory:
#         if ts.time % 10 == 0:
#             W.write(cluster)

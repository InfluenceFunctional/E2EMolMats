import MDAnalysis as mda
import os
import wandb
from cluster_figs import \
    (plot_rdf_series, plot_intermolecular_rdf_series, plot_cluster_stability)
from utils import dict2namespace

params = {
    'battery_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\battery_5/',
    'machine': 'local',  # or 'cluster'
    'show_figs': False,
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
        if run_dir != 'common':
            os.chdir(config.battery_path)
            os.chdir(run_dir)
            # do the analysis
            u = mda.Universe("system.data", "traj.dcd", format="LAMMPS")

            fig = plot_rdf_series(u)
            if config.show_figs:
                fig.show()
            wandb.log({'RDF Series': fig})

            fig = plot_intermolecular_rdf_series(u)
            if config.show_figs:
                fig.show()
            wandb.log({'Intermolecular RDF Series': fig})

            # todo atomwise intermolecular rdfs

            fig = plot_cluster_stability(u)
            if config.show_figs:
                fig.show()
            wandb.log({'Cluster Stability': fig})

    aa = 0



# cluster = u.select_atoms("all")
# with mda.Writer("traj.xyz", cluster.n_atoms) as W:
#     for ts in u.trajectory:
#         if ts.time % 10 == 0:
#             W.write(cluster)

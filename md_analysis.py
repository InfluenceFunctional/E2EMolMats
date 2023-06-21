import MDAnalysis as mda
import os
import wandb
from cluster_figs import \
    (plot_rdf_series, plot_cluster_stability)
from utils import dict2namespace

params = {
    'battery_path': r'C:\Users\mikem\crystals\clusters\cluster_structures\battery_4/',
}
config = dict2namespace(params)

with wandb.init(config=config, project="nicotinamide_clusters",
                entity="mkilgour", tags="reporting_test1"):

    wandb.run.name = wandb.config.machine + '_' + str(
        wandb.config.run_num)  # overwrite procedurally generated run name with our run name
    wandb.run.save()

    os.chdir(config.battery_path)
    dirs = os.listdir()
    for run_dir in dirs:  # loop over run directories in the battery
        if run_dir != 'common':
            # do the analysis
            u = mda.Universe("system.data", "traj.dcd", format="LAMMPS")

            fig = plot_rdf_series(u)
            wandb.log({'RDF Series': fig})

            fig = plot_cluster_stability(u)
            wandb.log({'Cluster Stability': fig})

    aa = 0

'''
Cluster analysis:
:> is the cluster stable (does it fall apart)
    --> need to assign molecules as 'inside' or 'outside' the cluster :: DONE
:> what is the intermolecular packing like
    --> compare intermolecular atomwise rdfs
'''

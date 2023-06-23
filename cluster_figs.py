import numpy as np
import plotly.graph_objects as go
from plotly.colors import n_colors
from sklearn.cluster import AgglomerativeClustering
from scipy.stats import mode
from MDAnalysis.analysis import rdf
import MDAnalysis as mda
from scipy.spatial import distance_matrix
import os


def plot_rdf_series(u):

    total_time = u.trajectory.totaltime
    n_steps = 10
    times = np.arange(0, total_time + 1, total_time // n_steps)

    rdf_analysis = rdf.InterRDF(u.atoms, u.atoms, range=(0.5, 6), nbins=200, verbose=False)
    n_frames = u.trajectory.n_frames
    rdf_step = n_frames // n_steps

    rdfs = []
    for step in range(0, n_frames, rdf_step):  # range(0, n_frames - 1, rdf_step):
        rdf_analysis.run(start=step, stop=step + 1, step=1, verbose=False)
        rdfs.append(rdf_analysis.results.rdf)
    rdfs = np.asarray(rdfs)

    colors = n_colors('rgb(250,50,5)', 'rgb(5,120,200)', len(rdfs), colortype='rgb')
    fig = go.Figure()
    for i in range(len(rdfs)):
        fig.add_trace(go.Scattergl(x=rdf_analysis.results.bins, y=rdfs[i], name=f'{times[i]:.0f} ps', marker=dict(color=colors[i])))

    fig.update_layout(xaxis_title='Range (A)', yaxis_title='Full RDF')

    return fig


def plot_intermolecular_rdf_series(u):
    n_molecules = 50
    randints = np.random.randint(0, len(u.residues), size=n_molecules)

    rdfs_list = []
    total_time = u.trajectory.totaltime
    n_steps = 10
    times = np.arange(0, total_time + 1, total_time // n_steps)
    for mol_ind in randints:
        mol = u.residues[mol_ind].atoms
        inter_mols = sum([u.residues[ind] for ind in range(len(u.residues)) if ind != mol_ind]).atoms
        rdf_analysis = rdf.InterRDF(mol, inter_mols, range=(0.5, 10), nbins=200, verbose=False)
        n_frames = u.trajectory.n_frames
        rdf_step = n_frames // n_steps

        rdfs = []
        for step in range(0, n_frames, rdf_step):  # range(0, n_frames - 1, rdf_step):
            rdf_analysis.run(start=step, stop=step + 1, step=1, verbose=False)
            rdfs.append(rdf_analysis.results.rdf)
        rdfs = np.asarray(rdfs)
        rdfs_list.append(rdfs)
    rdfs_list = np.asarray(rdfs_list)  # [molecules, time_steps, rdf_bins]
    combined_rdfs = rdfs_list.sum(0)  # average over molecules

    colors = n_colors('rgb(250,50,5)', 'rgb(5,120,200)', len(combined_rdfs), colortype='rgb')
    fig = go.Figure()
    for i in range(len(combined_rdfs)):
        fig.add_trace(go.Scattergl(x=rdf_analysis.results.bins, y=combined_rdfs[i], name=f'{times[i]:.0f} ps', marker=dict(color=colors[i])))

    fig.update_layout(xaxis_title='Range (A)', yaxis_title='Intermolecular RDF')

    return fig


def plot_cluster_stability(u: mda.Universe):
    mol_size = np.ptp(u.residues[0].atoms.positions)
    clusters_list_list = []
    majority_cluster_list = []
    majority_proportion_list = []
    radii = [4, 5, 6, 7, 8, 10, 12]

    ps_step = 100
    total_time = u.trajectory.totaltime
    times = np.arange(0, total_time + 1, ps_step)

    for cluster_threshold in radii:  # try different cutoffs
        '''identify molecules and assign to inside or outside cluster'''
        clusters_list = []
        majority_cluster = []
        for ts in u.trajectory:
            if ts.time % ps_step == 0:
                molecules = u.residues
                centroids = np.asarray([molecules[i].atoms.centroid() for i in range(len(molecules))])
                clustering = AgglomerativeClustering(linkage='single', metric='euclidean', distance_threshold=cluster_threshold, n_clusters=None).fit(centroids)
                clusters_list.append(clustering.labels_)
                modal, counts = mode(clustering.labels_, keepdims=False)
                majority_cluster.append(modal)
        clusters_list = np.asarray(clusters_list)
        majority_cluster = np.asarray(majority_cluster)
        clusters_list_list.append(clusters_list)
        majority_cluster_list.append(majority_cluster)

        majority_proportion = np.asarray([
            np.average(mols == majority) for mols, majority in zip(clusters_list, majority_cluster)
        ])

        majority_proportion_list.append(majority_proportion)
    #
    # clusters_list_list = np.asarray(clusters_list_list)
    # majority_cluster_list = np.asarray(majority_cluster_list)
    majority_proportion_list_list = np.asarray(majority_proportion_list)

    colors = n_colors('rgb(250,50,5)', 'rgb(5,120,200)', len(majority_proportion_list_list), colortype='rgb')

    fig = go.Figure()
    for i, radius in enumerate(radii):
        fig.add_trace(go.Scattergl(
            x=times, y=majority_proportion_list_list[i],
            name=f'Cutoff {radius:.1f} (A)', fill='tonexty', marker=dict(color=colors[i])
        ))
    fig.update_yaxes(range=[-0.05, 1.05])
    fig.update_layout(xaxis_title='Time (ps)', yaxis_title='Proportion in Majority Cluster')

    return fig


def plot_cluster_centroids_drift(u: mda.Universe):

    ps_step = 100
    total_time = u.trajectory.totaltime
    times = np.arange(0, total_time + 1, ps_step)

    distmat_list = []
    for ts in u.trajectory:
        if ts.time % ps_step == 0:
            molecules = u.residues
            centroids = np.asarray([molecules[i].atoms.centroid() for i in range(len(molecules))])
            distmat_list.append(distance_matrix(centroids, centroids))
    distmat_list = np.asarray(distmat_list)

    distmat_drift = np.zeros(len(distmat_list))
    for i in range(len(distmat_list)):
        distmat_drift[i] = np.abs(np.sum((distmat_list[i] - distmat_list[0]))) / distmat_list[0].sum()

    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=times, y=distmat_drift))
    fig.update_yaxes(range=[-0.05, 1.05])
    fig.update_layout(xaxis_title='Time (ps)', yaxis_title='Normed Intermolecular Centroids Drift')

    return fig


def process_thermo_data():

    f = open('screen.log', "r")
    text = f.read()
    lines = text.split('\n')
    f.close()

    skip = True
    results_dict = {'time step':[],
                    'temp':[],
                    'E_pair':[],
                    'E_mol':[],
                    'E_tot':[],
                    'Press':[]}
    for ind, line in enumerate(lines):
        if skip:
            if "Step" in line:
                skip = False
                #print(ind)
        else:
            if "Loop" in line:
                skip = True
                #print(ind)

            if not skip:
                split_line = line.split(' ')
                entries = [float(entry) for entry in split_line if entry != '']
                for ind2, key in enumerate(results_dict.keys()):
                    results_dict[key].append(entries[ind2])

    for key in results_dict.keys():
        results_dict[key] = np.asarray(results_dict[key])

    return results_dict

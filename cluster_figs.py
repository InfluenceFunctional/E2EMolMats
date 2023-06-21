import numpy as np
import plotly.graph_objects as go
from plotly.colors import n_colors
from sklearn.cluster import AgglomerativeClustering
from scipy.stats import mode
from MDAnalysis.analysis import rdf


def plot_rdf_series(u):
    rdf_analysis = rdf.InterRDF(u.atoms, u.atoms, range=(0.5, 6), nbins=200, verbose=True)
    n_frames = u.trajectory.n_frames
    rdf_step = n_frames // 10

    rdfs = []
    for step in range(0, 20):  # range(0, n_frames - 1, rdf_step):
        rdf_analysis.run(start=step, stop=step + 1, step=1, verbose=True)
        rdfs.append(rdf_analysis.results.rdf)
    rdfs = np.asarray(rdfs)

    colors = n_colors('rgb(250,50,5)', 'rgb(5,120,200)', len(rdfs), colortype='rgb')
    fig = go.Figure()
    for i in range(len(rdfs)):
        fig.add_trace(go.Scattergl(x=rdf_analysis.results.bins, y=rdfs[i], name=i, marker=dict(color=colors[i])))
    fig.show()
    return fig


def plot_cluster_stability(u: mda.Universe):
    mol_size = np.ptp(u.residues[0].atoms.positions)
    clusters_list_list = []
    majority_cluster_list = []
    majority_proportion_list = []
    radii = [4, 6, 8, 10, 12]
    for cluster_threshold in radii:  # try different cutoffs
        '''identify molecules and assign to inside or outside cluster'''
        clusters_list = []
        majority_cluster = []
        for ts in u.trajectory:
            if ts.time % 10 == 0:
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

    clusters_list_list = np.asarray(clusters_list_list)
    majority_cluster_list = np.asarray(majority_cluster_list)
    majority_proportion_list_list = np.asarray(majority_proportion_list)

    fig = go.Figure()
    for i, radius in enumerate(radii):
        fig.add_trace(go.Scattergl(x=np.arange(len(clusters_list)), y=majority_proportion_list_list[i], name=f'Cutoff {radius:.1f}'))
    fig.update_yaxes(range=[-0.05, 1.05])
    fig.show()
    return fig



# cluster = u.select_atoms("all")
# with mda.Writer("traj.xyz", cluster.n_atoms) as W:
#     for ts in u.trajectory:
#         if ts.time % 10 == 0:
#             W.write(cluster)
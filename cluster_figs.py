import numpy as np
import plotly.graph_objects as go
from plotly.colors import n_colors
from sklearn.cluster import AgglomerativeClustering
from scipy.stats import mode
from MDAnalysis.analysis import rdf
import MDAnalysis as mda
from scipy.spatial import distance_matrix
import os
from utils import compute_rdf_distance, compute_principal_axes_np, ff_names_dict
from plotly.subplots import make_subplots


def plot_rdf_series(u):

    total_time = u.trajectory.totaltime
    n_steps = 10
    times = np.arange(0, total_time + 1, total_time // n_steps)

    rdf_analysis = rdf.InterRDF(u.atoms, u.atoms, range=(0.5, 10), nbins=200, verbose=False, norm='rdf')
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

    return fig, rdf_analysis.results.bins, rdfs


def plot_intermolecular_rdf_series(u, atom_type1=None, atom_type2=None, n_frames_avg=1):
    n_molecules = min((100, len(u.residues)))
    randints = np.random.choice(len(u.residues), size=n_molecules, replace=False)

    rdfs_list = []
    total_time = u.trajectory.totaltime
    n_steps = 10
    times = np.arange(0, total_time + 1, total_time // n_steps)
    for mol_ind in randints:
        mol = u.residues[mol_ind].atoms
        inter_mols = sum([u.residues[ind] for ind in range(len(u.residues)) if ind != mol_ind]).atoms

        if atom_type1 is not None:
            mol = mol.select_atoms("name " + atom_type1)
        if atom_type2 is not None:
            inter_mols = inter_mols.select_atoms("name " + atom_type2)

        rdf_analysis = rdf.InterRDF(mol, inter_mols, range=(0.5, 10), nbins=200, verbose=False, norm='rdf')
        n_frames = u.trajectory.n_frames
        rdf_step = n_frames // n_steps

        rdfs = []
        for step in range(0, n_frames, rdf_step):  # range(0, n_frames - 1, rdf_step):
            rdf_analysis.run(start=step, stop=step + n_frames_avg, step=1, verbose=False)
            rdfs.append(rdf_analysis.results.rdf)
        rdfs = np.asarray(rdfs)
        rdfs_list.append(rdfs)
    rdfs_list = np.asarray(rdfs_list)  # [molecules, time_steps, rdf_bins]
    combined_rdfs = rdfs_list.mean(0)  # average over molecules

    colors = n_colors('rgb(250,50,5)', 'rgb(5,120,200)', len(combined_rdfs), colortype='rgb')
    fig = go.Figure()
    for i in range(len(combined_rdfs)):
        fig.add_trace(go.Scattergl(x=rdf_analysis.results.bins, y=combined_rdfs[i], name=f'{times[i]:.0f} ps', marker=dict(color=colors[i])))

    fig.update_layout(xaxis_title='Range (A)', yaxis_title='Intermolecular RDF')

    return fig, rdf_analysis.results.bins, combined_rdfs


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

    return fig, times, majority_proportion_list_list


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
        distmat_drift[i] = np.sum(np.abs((distmat_list[i] - distmat_list[0]))) / distmat_list[0].sum()

    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=times, y=distmat_drift))
    fig.update_yaxes(range=[-0.05, max(1.05, max(distmat_drift))])
    fig.update_layout(xaxis_title='Time (ps)', yaxis_title='Normed Intermolecular Centroids Drift')

    return fig, times, distmat_drift


def process_thermo_data():

    f = open('screen.log', "r")
    text = f.read()
    lines = text.split('\n')
    f.close()

    skip = True
    results_dict = {'time step': [],
                    'temp': [],
                    'E_pair': [],
                    'E_mol': [],
                    'E_tot': [],
                    'Press': []}
    for ind, line in enumerate(lines):
        if skip:
            if "Step" in line:
                skip = False
                # print(ind)
        else:
            if "Loop" in line:
                skip = True
                # print(ind)

            if not skip:
                split_line = line.split(' ')
                entries = [float(entry) for entry in split_line if entry != '']
                for ind2, key in enumerate(results_dict.keys()):
                    results_dict[key].append(entries[ind2])

    for key in results_dict.keys():
        results_dict[key] = np.asarray(results_dict[key])

    return results_dict


def plot_atomwise_rdf_drift(u, atomwise_rdfs, bins):
    t0_atomwise_rdfs = np.asarray([
        rdf[0] for rdf in atomwise_rdfs
    ])
    num_traj_points = len(atomwise_rdfs[0])
    trajectory_atomwise_rdfs = [
        np.asarray([
            rdf[i] for rdf in atomwise_rdfs
        ]) for i in range(1, num_traj_points)
    ]
    rdf_drift = np.zeros(num_traj_points)
    for i in range(1, num_traj_points):
        rdf_drift[i] = compute_rdf_distance(t0_atomwise_rdfs, trajectory_atomwise_rdfs[i - 1], bins)

    ps_step = 100
    total_time = u.trajectory.totaltime
    times = np.arange(0, total_time + 1, ps_step)

    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=times, y=rdf_drift))
    fig.update_yaxes(range=[-0.05, max(1.05, max(rdf_drift))])
    fig.update_layout(xaxis_title='Time (ps)', yaxis_title='Intermolecular Atomwise RDF Drift')

    return fig, times, rdf_drift


def plot_atomwise_rdf_ref_dist(u, atomwise_rdfs, ref_atomwise_rdfs, bins):
    num_traj_points = len(atomwise_rdfs[0])
    trajectory_atomwise_rdfs = [
        np.asarray([
            rdf[i] for rdf in atomwise_rdfs
        ]) for i in range(num_traj_points)
    ]

    num_ref_points = len(ref_atomwise_rdfs[0]) - 1 # give the reference time to equilibrade
    reference_atomwise_rdfs = [
        np.asarray([
            rdf[i] for rdf in ref_atomwise_rdfs
        ]) for i in range(1, num_ref_points)
    ]

    rdf_drift = np.zeros((num_traj_points, num_ref_points))
    for i in range(num_traj_points):
        for j in range(num_ref_points - 1):
            rdf_drift[i, j] = compute_rdf_distance(reference_atomwise_rdfs[j], trajectory_atomwise_rdfs[i], bins)

    mean_rdf_drift = rdf_drift.mean(1)  # average over reference trajectory

    ps_step = 100
    total_time = u.trajectory.totaltime
    times = np.arange(0, total_time + 1, ps_step)

    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=times, y=mean_rdf_drift))
    fig.update_yaxes(range=[-0.05, max(1.05, max(mean_rdf_drift))])
    fig.update_layout(xaxis_title='Time (ps)', yaxis_title='Intermolecular Atomwise RDF Drift')

    return fig, times, mean_rdf_drift


def plot_alignment_fingerprint(u):
    ps_step = 100
    total_time = u.trajectory.totaltime
    times = np.arange(0, total_time + 1, ps_step)

    Ip_overlaps_list = []
    for ts in u.trajectory:
        if ts.time % ps_step == 0:
            molecules = u.residues
            coords = np.asarray([molecules[i].atoms.positions for i in range(len(molecules))])
            Ip_list = []
            for j in range(len(molecules)):
                Ip, _, _ = compute_principal_axes_np(coords[j])
                Ip_list.append(Ip)
            Ip_list = np.stack(Ip_list)

            # source mol, target mol, source Ip, target Ip
            Ip_overlaps = np.zeros((len(molecules), len(molecules), 3, 3))
            for j in range(len(molecules)):
                for k in range(3):
                    Ip_overlaps[j, :, k, :] = Ip_list[j, k].dot(np.transpose(Ip_list, axes=[0, 2, 1]))

            Ip_overlaps_list.append(Ip_overlaps)

    Ip_overlaps_list = np.stack(Ip_overlaps_list)

    Ip_overlaps_drift = np.zeros(len(Ip_overlaps_list))
    for i in range(len(Ip_overlaps_list)):
        Ip_overlaps_drift[i] = np.sum(np.abs((Ip_overlaps_list[i] - Ip_overlaps_list[0]))) / np.prod(list(Ip_overlaps_list[0].shape))  # norm by maximum possible values

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=times, y=Ip_overlaps_drift))
    fig.update_yaxes(range=[-0.05, max(1.05, max(Ip_overlaps_drift))])
    fig.update_layout(xaxis_title='Time (ps)', yaxis_title='Normed Molecule Principal Axes Overlap Drift')

    return fig, Ip_overlaps_drift, times


def plot_thermodynamic_data(thermo_results_dict):
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
    return fig


def trajectory_rdf_analysis(u):
    ""
    '''full rdf'''
    atom_types = u.atoms.types
    atom_names = np.asarray([ff_names_dict[atype] for atype in atom_types])
    u.add_TopologyAttr('name', atom_names)
    fig, bins, full_rdf = plot_rdf_series(u)
    # if config.show_figs:
    #     fig.show(renderer="browser")
    # wandb.log({'RDF Series': fig})
    # fig.write_image('RDF_series.png')

    '''intermolecular df'''
    fig, bins, intermolecular_rdf = plot_intermolecular_rdf_series(u)
    # if config.show_figs:
    #     fig.show(renderer="browser")
    # wandb.log({'Intermolecular RDF Series': fig})
    # fig.write_image('intermolecular_RDF_series.png')

    '''atom-type-wise intermolecular rdf'''
    atomwise_rdfs = []
    for atom_type in ff_names_dict.values():
        if 'h' not in atom_type:
            fig, bins, atom_type_rdf = plot_intermolecular_rdf_series(u, atom_type, atom_type, n_frames_avg=1)
            atomwise_rdfs.append(atom_type_rdf)

    #         if config.show_figs:
    #             fig.show(renderer="browser")
    #         wandb.log({atom_type + ' Intermolecular RDF Series': fig})
    #         fig.write_image(atom_type + ' intermolecular_RDF_series.png')
    #
    return full_rdf, intermolecular_rdf, atomwise_rdfs, bins

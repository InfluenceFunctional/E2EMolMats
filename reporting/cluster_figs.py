import numpy as np
import plotly.graph_objects as go
from plotly.colors import n_colors
from sklearn.cluster import AgglomerativeClustering
from scipy.stats import mode
from MDAnalysis.analysis import rdf
import MDAnalysis as mda
from scipy.spatial import distance_matrix
import os
from utils import compute_rdf_distance, compute_principal_axes_np, ff_names_dict, tile_universe
from plotly.subplots import make_subplots
from scipy.spatial.distance import cdist, pdist
from analysis.cluster_rdf import crystal_rdf


def plot_rdf_series(u, nbins, rrange, core_inds, rdf_norm='rdf', n_frames_avg=1):
    total_time = u.trajectory.totaltime
    n_steps = 10
    times = np.arange(0, total_time + 1, total_time // n_steps)

    rdf_analysis = rdf.InterRDF(u.residues[core_inds].atoms, u.atoms, range=rrange, nbins=nbins, verbose=False, norm=rdf_norm)
    n_frames = u.trajectory.n_frames
    rdf_step = n_frames // n_steps

    rdfs = []
    for step in range(0, n_frames, rdf_step):  # range(0, n_frames - 1, rdf_step):
        rdf_analysis.run(start=step, stop=step + n_frames_avg, step=1, verbose=False)
        rdfs.append(rdf_analysis.results.rdf)
    rdfs = np.asarray(rdfs)

    colors = n_colors('rgb(250,50,5)', 'rgb(5,120,200)', len(rdfs), colortype='rgb')
    fig = go.Figure()
    for i in range(len(rdfs)):
        fig.add_trace(go.Scattergl(x=rdf_analysis.results.bins, y=rdfs[i], name=f'{times[i]:.0f} ps', marker=dict(color=colors[i])))

    fig.update_layout(xaxis_title='Range (A)', yaxis_title='Full RDF')

    return fig, rdf_analysis.results.bins, rdfs


def plot_intermolecular_rdf_series(u, nbins, rrange, core_inds, atom_type1=None, atom_type2=None, rdf_norm='rdf', n_frames_avg=1):
    n_molecules = min((100, len(core_inds)))
    randints = np.random.choice(core_inds, size=n_molecules, replace=False)

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

        rdf_analysis = rdf.InterRDF(mol, inter_mols, range=rrange, nbins=nbins, verbose=False, norm=rdf_norm)
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
    rdf_drift = np.zeros((len(atomwise_rdfs), len(ref_atomwise_rdfs)))
    for i in range(len(atomwise_rdfs)):
        for j in range(len(ref_atomwise_rdfs)):
            rdf_drift[i, j] = compute_rdf_distance(atomwise_rdfs[i], ref_atomwise_rdfs[j], bins, envelope='tanh')

    mean_rdf_drift = rdf_drift.mean(1)  # average over reference trajectory

    return mean_rdf_drift


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


def trajectory_rdf_analysis(u, n_mols_to_sample=10, nbins=200, rrange=[0, 10], core_cutoff=None, tiling=None):
    ""
    '''update atom types'''
    atom_types = u.atoms.types
    atom_names = np.asarray([ff_names_dict[atype] for atype in atom_types])
    u.add_TopologyAttr('name', atom_names)
    mol_num_atoms = u.residues[0].atoms.n_atoms
    total_time = u.trajectory.totaltime
    n_frames = u.trajectory.n_frames
    time_step = u.trajectory.dt
    print_steps = 10

    print_frames = np.arange(n_frames // print_steps, n_frames + 1, step=n_frames // print_steps)

    frames = []
    rdfs = []
    for ts in u.trajectory:
        if ts.frame in print_frames:
            frames.append(ts.frame)
            if tiling is not None:
                new_u = tile_universe(u, tiling)
            else:
                new_u = u

            all_inds = np.arange(len(new_u.atoms))
            if core_cutoff is not None:
                molecules = new_u.residues
                centroids = np.asarray([molecules[i].atoms.centroid() for i in range(len(molecules))])
                global_centroid = centroids.mean(0)
                dists = cdist(global_centroid[None, :], centroids)[0]
                core_mols = np.argwhere(dists < core_cutoff)[:, 0]
            else:
                core_mols = np.random.choice(len(u.residues), size=n_mols_to_sample)

            atomwise_rdfs = np.zeros(((mol_num_atoms ** 2 - mol_num_atoms) // 2 + mol_num_atoms, nbins))
            n_samples = min(n_mols_to_sample, len(core_mols))
            for mm in range(n_samples):
                core_inds = new_u.residues[core_mols[mm]].atoms.indices  # one inside molecule at a time
                outer_inds = np.concatenate([new_u.residues[ii].atoms.indices for ii in range(len(new_u.residues)) if ii != core_mols[mm]])

                atomwise_rdfs_i, bins, rdfs_dict = crystal_rdf(
                    positions=new_u.atoms.positions,
                    symbols=new_u.atoms.types.astype(np.int64),
                    mol_num_atoms=mol_num_atoms,
                    in_inds=core_inds,
                    out_inds=outer_inds,
                    rrange=rrange,
                    bins=nbins,
                    raw_density=True,
                    elementwise=False,
                    atomwise=True,
                )
                atomwise_rdfs += atomwise_rdfs_i[0] / n_samples  # manual average

            rdfs.append(atomwise_rdfs)

    rdfs = np.stack(rdfs)

    return rdfs, bins, np.asarray(frames) * time_step


def old_trajectory_rdf_analysis(u, nbins=200, rrange=[0, 10], core_cutoff=None, rdf_norm='rdf', n_frames_avg=1):
    ""
    if core_cutoff is not None:  # compute RDF for 'core' molecules
        molecules = u.residues
        centroids = np.asarray([molecules[i].atoms.centroid() for i in range(len(molecules))])
        global_centroid = centroids.mean(0)
        dists = cdist(global_centroid[None, :], centroids)[0]
        core_inds = np.argwhere(dists < core_cutoff)[:, 0]
    else:
        core_inds = np.arange(len(u.residues))

    '''full rdf'''
    atom_types = u.atoms.types
    atom_names = np.asarray([ff_names_dict[atype] for atype in atom_types])
    u.add_TopologyAttr('name', atom_names)
    fig, nbins, full_rdf = plot_rdf_series(u, nbins, rrange, core_inds, rdf_norm=rdf_norm, n_frames_avg=n_frames_avg)

    '''intermolecular df'''
    fig, nbins, intermolecular_rdf = plot_intermolecular_rdf_series(u, nbins, rrange, core_inds, rdf_norm=rdf_norm, n_frames_avg=n_frames_avg)

    '''atom-type-wise intermolecular rdf'''
    atomwise_rdfs = []
    for atom_type in ff_names_dict.values():
        if 'h' not in atom_type:
            fig, nbins, atom_type_rdf = plot_intermolecular_rdf_series(u, nbins, rrange, core_inds, atom_type, atom_type, rdf_norm=rdf_norm, n_frames_avg=n_frames_avg)
            atomwise_rdfs.append(atom_type_rdf)

    return full_rdf, intermolecular_rdf, atomwise_rdfs, nbins

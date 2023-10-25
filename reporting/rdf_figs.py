import numpy as np
from scipy.spatial.distance import cdist

from analysis.cluster_rdf import crystal_rdf
from reporting.cluster_figs import plot_rdf_series, plot_intermolecular_rdf_series
from utils import ff_names_dict, tile_universe


def trajectory_rdf_analysis(u, n_mols_to_sample=10, nbins=200, rrange=[0, 10], core_cutoff=None, tiling=None, print_steps=10):
    ""
    '''update atom types'''
    atom_types = u.atoms.types
    atom_names = np.asarray([ff_names_dict[atype] for atype in atom_types])
    u.add_TopologyAttr('name', atom_names)
    mol_num_atoms = u.residues[0].atoms.n_atoms
    total_time = u.trajectory.totaltime
    n_frames = u.trajectory.n_frames
    time_step = u.trajectory.dt

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

            mol_num_elements = 9
            atomwise_rdfs = np.zeros(((mol_num_elements ** 2 - mol_num_elements) // 2 + mol_num_elements, nbins))  # np.zeros(((mol_num_atoms ** 2 - mol_num_atoms) // 2 + mol_num_atoms, nbins))
            n_samples = min(n_mols_to_sample, len(core_mols))
            for mm in range(n_samples):
                core_inds = new_u.residues[core_mols[mm]].atoms.indices  # one inside molecule at a time
                outer_inds = np.concatenate([new_u.residues[ii].atoms.indices for ii in range(len(new_u.residues)) if ii != core_mols[mm]])  # limit to 15 atoms - kill trailing benzamide proton

                atomwise_rdfs_i, bins, rdfs_dict = crystal_rdf(
                    positions=new_u.atoms.positions,
                    symbols=new_u.atoms.types.astype(np.int64),
                    mol_num_atoms=mol_num_atoms,
                    in_inds=core_inds,
                    out_inds=outer_inds,
                    rrange=rrange,
                    bins=nbins,
                    raw_density=True,
                    elementwise=True,
                    atomwise=False,
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

from argparse import Namespace
import os
import numpy as np
import MDAnalysis as mda
from scipy.spatial.distance import cdist, pdist

names_dict = {'1': 'H',  # rename for xyz export
              '8': 'H',
              '2': 'H',
              '6': 'N',
              '7': 'N',
              '4': 'C',
              '5': 'C',
              '3': 'O',
              }

ff_names_dict = {'1': 'ha',  # detailed atom types for analysis
                 '8': 'h4',
                 '2': 'hn',
                 '6': 'n',
                 '7': 'nb',
                 '4': 'c',
                 '5': 'ca',
                 '3': 'o',
                 }


def dict2namespace(data_dict: dict):
    """
    Recursively converts a dictionary and its internal dictionaries into an
    argparse.Namespace

    Parameters
    ----------
    data_dict : dict
        The input dictionary

    Return
    ------
    data_namespace : argparse.Namespace
        The output namespace
    """
    for k, v in data_dict.items():
        if isinstance(v, dict):
            data_dict[k] = dict2namespace(v)
        else:
            pass
    data_namespace = Namespace(**data_dict)

    return data_namespace


def rotvec2sph(rotvec):
    """
    transform rotation vector with axis rotvec/norm(rotvec) and angle ||rotvec||
    to spherical coordinates theta, phi and r ||rotvec||
    """
    r = np.linalg.norm(rotvec, axis=-1)
    if rotvec.ndim == 1:
        rotvec = rotvec[None, :]
        r = np.asarray(r)[None]

    unit_vector = rotvec / r[:, None]

    # convert unit vector to angles
    theta = np.arctan2(np.sqrt(unit_vector[:, 0] ** 2 + unit_vector[:, 1] ** 2), unit_vector[:, 2])
    phi = np.arctan2(unit_vector[:, 1], unit_vector[:, 0])
    if rotvec.ndim == 1:
        return np.concatenate((theta, phi, r), axis=-1)  # polar, azimuthal, applied rotation
    else:
        return np.concatenate((theta[:, None], phi[:, None], r[:, None]), axis=-1)  # polar, azimuthal, applied rotation


def compute_Ip_alignment(u, Ip_overlap_trajectory):
    cutoffs = [0.75, 0.8, 0.85, 0.9, 0.95]
    Ip_alignment_trajectory = np.zeros((len(Ip_overlap_trajectory), len(cutoffs)))
    for ic in range(len(cutoffs)):
        alignment_cutoff = cutoffs[ic]
        for it in range(len(Ip_overlap_trajectory)):
            alignments = np.zeros((3, len(u.residues), len(u.residues)))
            for ip in range(3):
                alignments[ip] = Ip_overlap_trajectory[it][:, :, ip, ip] > alignment_cutoff
            tot_alignment = np.prod(alignments, axis=0)
            Ip_alignment_trajectory[it, ic] = tot_alignment.mean()

    return Ip_alignment_trajectory


def compute_Ip_molwise_alignment(u, Ip_overlap_trajectory, print_frames):
    cutoffs = [0.25, 0.5, 0.75, 0.9, 0.95]

    global_Ip_alignment_trajectory = np.zeros((len(Ip_overlap_trajectory), u.residues.n_residues, len(cutoffs)))
    local_Ip_alignment_trajectory = np.zeros((len(Ip_overlap_trajectory), u.residues.n_residues))

    it = -1
    for itt, ts in enumerate(u.trajectory):
        if itt in print_frames:
            it += 1
            mol_centroids = np.asarray([u.residues[ii].atoms.centroid() for ii in range(len(u.residues))])
            distmat = cdist(mol_centroids, mol_centroids)
            envelope = -np.tanh((distmat - 10) / 2) / 2 + .5  # a function which prioritizes mols which are nearby to each other

            for icc, cutoff in enumerate(cutoffs):
                alignments = np.zeros((3, len(u.residues), len(u.residues)))

                for ip in range(3):
                    alignments[ip] = Ip_overlap_trajectory[it, :, :, ip, ip]

                tot_alignment = np.mean(alignments > cutoff, axis=0)
                global_Ip_alignment_trajectory[it, :, icc] = (tot_alignment).mean(1)

            local_alignments = np.zeros((3, len(u.residues), len(u.residues)))
            for ip in range(3):
                local_alignments[ip] = np.abs(Ip_overlap_trajectory[it, :, :, ip, ip]) * envelope

            local_Ip_alignment_trajectory[it, :] = local_alignments.mean(0).sum(1)

    return global_Ip_alignment_trajectory, local_Ip_alignment_trajectory


def compute_rdf_distance(rdf1, rdf2, rr, envelope=None):
    """
    compute a distance metric between two radial distribution functions with shapes
    [num_sub_rdfs, num_bins] where sub_rdfs are e.g., particular interatomic RDFS within a certain sample (elementwise or atomwise modes)
    rr is the bin edges used for both rdfs

    option for input to be torch tensors or numpy arrays, but has to be the same either way
    computation is done in torch
    range rr can be independently either np.array or torch.tensor
    will return same format as given
    """
    if envelope is None:
        tempering_func = np.ones(rdf1.shape)
    elif envelope == 'tanh':
        x = np.linspace(-10, 2, rdf1.shape[-1])
        tempering_func = (np.tanh(-x) / 2) + 0.5
        tempering_func = tempering_func[None, :]

    smoothed_rdf1 = rdf1 * tempering_func
    smoothed_rdf2 = rdf2 * tempering_func

    emd = earth_movers_distance_np(smoothed_rdf1, smoothed_rdf2)

    range_normed_emd = emd * (rr[1] - rr[0]) ** 2  # rescale the distance from units of bins to the real physical range

    distance = range_normed_emd.mean()  # rescale by respective densities?

    assert np.sum(np.isnan(distance)) == 0

    return distance


def earth_movers_distance_np(d1: np.ndarray, d2: np.ndarray):
    """
    earth mover's distance between two PDFs
    not normalized or aggregated
    """
    return np.sum(np.abs(np.cumsum(d1, axis=-1) - np.cumsum(d2, axis=-1)), axis=-1)


def compute_principal_axes_np(coords):
    """
    compute the principal axes for a given set of particle coordinates, ignoring particle mass
    use our overlap rules to ensure a fixed direction for all axes under almost all circumstances
    """  # todo harmonize with torch version - currently disagrees ~0.5% of the time
    points = coords - coords.mean(0)

    x, y, z = points.T
    Ixx = np.sum((y ** 2 + z ** 2))
    Iyy = np.sum((x ** 2 + z ** 2))
    Izz = np.sum((x ** 2 + y ** 2))
    Ixy = -np.sum(x * y)
    Iyz = -np.sum(y * z)
    Ixz = -np.sum(x * z)
    I = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])  # inertial tensor
    Ipm, Ip = np.linalg.eig(I)  # principal inertial tensor
    Ipm, Ip = np.real(Ipm), np.real(Ip)
    sort_inds = np.argsort(Ipm)
    Ipm = Ipm[sort_inds]
    Ip = Ip.T[sort_inds]  # want eigenvectors to be sorted row-wise (rather than column-wise)

    # cardinal direction is vector from CoM to the farthest atom
    dists = np.linalg.norm(points, axis=1)
    max_ind = np.argmax(dists)
    max_equivs = np.argwhere(np.round(dists, 8) == np.round(dists[max_ind], 8))[:, 0]  # if there are multiple equidistant atoms - pick the one with the lowest index
    max_ind = int(np.amin(max_equivs))
    direction = points[max_ind]
    direction = np.divide(direction, np.linalg.norm(direction))
    overlaps = Ip.dot(direction)  # check if the principal components point towards or away from the CoG
    signs = np.sign(overlaps)  # returns zero for zero overlap, but we want it to default to +1 in this case
    signs[signs == 0] = 1

    Ip = (Ip.T * signs).T  # if the vectors have negative overlap, flip the direction
    if np.any(np.abs(overlaps) < 1e-3):  # if any overlaps are vanishing, determine the direction via the RHR (if two overlaps are vanishing, this will not work)
        # align the 'good' vectors
        fix_ind = np.argmin(np.abs(overlaps))  # vector with vanishing overlap
        if compute_Ip_handedness(Ip) < 0:  # make sure result is right handed
            Ip[fix_ind] = -Ip[fix_ind]

    return Ip, Ipm, I


def compute_Ip_handedness(Ip):
    """
    determine the right or left handedness from the cross products of principal inertial axes
    """
    if Ip.ndim == 2:
        return np.sign(np.dot(Ip[0], np.cross(Ip[1], Ip[2])).sum())
    elif Ip.ndim == 3:
        return np.sign(np.dot(Ip[:, 0], np.cross(Ip[:, 1], Ip[:, 2], axis=1).T).sum(1))


def cell_vol(v, a, units='natural'):
    ''' Calculate cos and sin of cell angles '''
    if units == 'degrees':
        a = a * np.pi / 180
    else:
        pass
    cos_a = np.cos(a)

    ''' Calculate volume of the unit cell '''
    val = 1.0 - cos_a[0] ** 2 - cos_a[1] ** 2 - cos_a[2] ** 2 + 2.0 * cos_a[0] * cos_a[1] * cos_a[2]
    vol = v[0] * v[1] * v[2] * np.sqrt(np.abs(val))  # technically a signed quanitity

    return vol


def rewrite_trajectory(u: mda.Universe, run_dir: str, extra_atomwise_values=None):

    #if not os.path.exists(f"{run_dir}_traj.xyz"):
    # if extra_atomwise_values is not None:
    #     u.add_TopologyAttr('tempfactors')
    atom_types = u.atoms.types
    atom_names = np.asarray([names_dict[atype] for atype in atom_types])
    u.add_TopologyAttr('name', atom_names)
    cluster = u.select_atoms("all")
    with mda.Writer(f"{run_dir}_traj.xyz", cluster.n_atoms) as W:
        for ts in u.trajectory:
            # if extra_atomwise_values is not None:
            #     u.atoms.tempfactors = extra_atomwise_values[ts.frame]
            W.write(cluster)

        if extra_atomwise_values is not None:
            newFile = open(f"{run_dir}_traj2.xyz", 'w')
            with open(f"{run_dir}_traj.xyz") as f:
                newText = f.readlines()
                counter = 0
                frame_num = None
                for ind, line in enumerate(newText):
                    if 'frame' in line:
                        newFile.write(line)
                        frame_num = int(line.split()[1])
                        counter = 0
                    elif len(line.split()) == 1:
                        newFile.write(line)
                    elif line == '\n':
                        newFile.write(line)
                    else:
                        coord_line = line.replace('\n', f'  {extra_atomwise_values[frame_num, counter]:.2f}\n')
                        counter += 1
                        newFile.write(coord_line)


def tile_universe(universe, tiling):
    n_x, n_y, n_z = tiling
    box = universe.dimensions[:3]
    copied = []
    for x in range(n_x):
        for y in range(n_y):
            for z in range(n_z):
                u_ = universe.copy()
                move_by = box * (x, y, z)
                u_.atoms.translate(move_by)
                copied.append(u_.atoms)

    new_universe = mda.Merge(*copied)
    new_box = box * (n_x, n_y, n_z)
    new_universe.dimensions = list(new_box) + [90] * 3
    return new_universe

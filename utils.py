from argparse import Namespace
from generate_cluster_structures import generate_structure
from ovito.io import import_file, export_file
from template_scripts.initial_setup_for_ovito import initial_setup
from template_scripts.original_templify_to_runnable import templify_to_runnable
import ovito
from distutils.dir_util import copy_tree
import os
import numpy as np

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


def create_xyz_and_run_lammps(head_dir, run_num, crystals_path, cluster_size,
                              cluster_type="supercell", structure_identifier="NICOAM13",
                              max_sphere_radius=None,
                              defect_rate=0, scramble_rate=0, gap_rate=0,
                              seed=1, min_inter_cluster_distance=200,
                              temperature=300, run_time=int(1e6),
                              print_steps=100, box_type='s',
                              integrator='langevin', damping: str = str(100.0)):
    """
    :param head_dir:
    :param run_num:
    :param crystals_path:
    :param cluster_size:
    :param cluster_type:
    :param structure_identifier:
    :param max_sphere_radius:
    :param defect_rate:
    :param scramble_rate:
    :param gap_rate:
    :param seed:
    :param min_inter_cluster_distance:
    :param temperature:
    :param run_time:
    :param print_steps:
    :param box_type:
    :param integrator:
    :param damping:
    :return:
    """

    '''make new workdir'''
    workdir = head_dir + '/' + str(run_num)
    if workdir is not None:
        if not os.path.exists(workdir):
            os.mkdir(workdir)
        os.chdir(workdir)
        '''copy in common elements'''
        copy_tree('../common', './')
    else:
        os.chdir(workdir)

    '''set temperature, run time, and print step in lmp file'''
    with open("run_MD.lmp") as f:
        newText = f.read().replace('_TEMP', str(temperature))
        newText = newText.replace('_RUNTIME', str(run_time))
        newText = newText.replace('_PRINTSTEPS', str(print_steps))
        newText = newText.replace('_SEED', str(seed))
        newText = newText.replace('_BOUND', str(box_type))
        newText = newText.replace('_DAMP', damping)
        if integrator == 'langevin':
            newText = newText.replace('#_LANGEVIN', '')
        elif integrator == 'nosehoover':
            newText = newText.replace('#_NOSE', '')
        elif integrator == 'npt':
            newText = newText.replace('#_NPT', '')

    with open("run_MD.lmp", "w") as f:
        f.write(newText)

    '''generate cluster structure'''
    xyz_filename = generate_structure(
        workdir, crystals_path, structure_identifier,
        cluster_type, max_sphere_radius,
        cluster_size, defect_rate, scramble_rate,
        gap_rate, seed, min_inter_cluster_distance,
        periodic_structure=box_type == 'p')

    '''convert from .xyz to lammps datafile'''
    pipeline = import_file(xyz_filename)
    export_file(pipeline, '1.data', 'lammps/data', atom_style='full')

    '''prep for ovito bonds'''
    initial_setup(workdir, '1.data', '2.data')

    '''add bonds via ovito'''
    ovito.scene.load("nicotinamide_bond_session.ovito")
    pipeline = ovito.scene.pipelines[0]
    pipeline.source.load('2.data')
    export_file(pipeline, '3.data', 'lammps/data', atom_style='full')

    '''ltemplify'''
    os.system('ltemplify.py 3.data > 4.lt')  # .py on ltemplify required on cluster not windows

    '''make runnable'''
    templify_to_runnable(workdir, "4.lt", "3.data", "5.lt")

    '''run moltemplate and cleanup'''
    os.system("moltemplate.sh system.lt")
    os.system("cleanup_moltemplate.sh")
    #
    # '''optionally - directly run MD'''
    os.system("sbatch sub_job.slurm")


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

    range_normed_emd = emd * (rr[1] - rr[0])  # rescale the distance from units of bins to the real physical range

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

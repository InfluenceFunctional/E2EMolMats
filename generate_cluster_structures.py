"""
script for converting unit cells to finite size crystal structures
options for cluster size and shape
options for introduction of defects
"""

import os
import sys

import numpy as np
from ase import Atoms
#from ase.visualize import view
from ase import io
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from utils import dict2namespace

def get_crystal_properties(structure_identifier):
    # get original structure
    if structure_identifier == "nicotinamide/NICOAM07":  # form 5, Tm 383.5, epsilon
        atoms_in_molecule = 15
        space_group = "P21"
        z_value = 2
        z_prime = 1
    elif structure_identifier == "nicotinamide/NICOAM08":  # form 7, Tm 381, eta
        atoms_in_molecule = 15
        space_group = "P-1"
        z_value = 4
        z_prime = 2
    elif structure_identifier == "nicotinamide/NICOAM09":  # form 8, Tm 377.5, theta
        atoms_in_molecule = 15
        space_group = "P21"
        z_value = 40
        z_prime = 20
    elif structure_identifier == "nicotinamide/NICOAM13":  # form 1, Tm 402K, alpha
        atoms_in_molecule = 15
        space_group = "P21/c"
        z_value = 4
        z_prime = 1
    elif structure_identifier == "nicotinamide/NICOAM14":  # form 2, Tm 390K, beta
        atoms_in_molecule = 15
        space_group = "P2/n"
        z_value = 16
        z_prime = 4
    elif structure_identifier == "nicotinamide/NICOAM15":  # form 3, Tm 388K, gamma
        atoms_in_molecule = 15
        space_group = "P21/c"
        z_value = 16
        z_prime = 4
    elif structure_identifier == "nicotinamide/NICOAM16":  # form 4, Tm 387K, delta
        atoms_in_molecule = 15
        space_group = "P21/c"
        z_value = 8
        z_prime = 2
    elif structure_identifier == "nicotinamide/NICOAM17":  # form 9, Tm 376K, iota
        atoms_in_molecule = 15
        space_group = "P21/c"
        z_value = 4
        z_prime = 1
    elif structure_identifier == "nicotinamide/NICOAM18":  # form 6, Tm 382.5K, zeta
        atoms_in_molecule = 15
        space_group = "P-1"
        z_value = 4
        z_prime = 2

    elif structure_identifier == "acridine/Form2":
        atoms_in_molecule = 23
        space_group = "P21/n"
        z_value = 4
        z_prime = 1
    elif structure_identifier == "acridine/Form3":
        atoms_in_molecule = 23
        space_group = "P21/c"
        z_value = 8
        z_prime = 2
    elif structure_identifier == "acridine/Form4":
        atoms_in_molecule = 23
        space_group = "P212121"
        z_value = 12
        z_prime = 3
    elif structure_identifier == "acridine/Form6":
        atoms_in_molecule = 23
        space_group = "Cc"
        z_value = 8
        z_prime = 2
    elif structure_identifier == "acridine/Form7":
        atoms_in_molecule = 23
        space_group = "P21/n"
        z_value = 8
        z_prime = 2
    elif structure_identifier == "acridine/Form8":
        atoms_in_molecule = 23
        space_group = "Cc"
        z_value = 8
        z_prime = 2
    elif structure_identifier == "acridine/Form9":
        atoms_in_molecule = 23
        space_group = "P21/n"
        z_value = 4
        z_prime = 1
    else:
        print("no such structure!")
        sys.exit()

    return atoms_in_molecule, space_group, z_value, z_prime
    # '''
    # II	04
    # III	07
    # IV	08
    # VI	05
    # VII	06
    # VIII	23bek9s	from Alex Shtukenberg
    # IX	12
    # '''


def generate_structure(workdir, crystals_path, structure_identifier,
                       cluster_type, max_sphere_radius, cluster_size,
                       defect_rate, scramble_rate, gap_rate, seed,
                       min_inter_cluster_distance, min_lattice_length,
                       periodic_structure=False, prep_crystal_in_melt=False):
    np.random.seed(seed=seed)

    # move to working directory
    if workdir is not None:
        os.chdir(workdir)

    atoms_in_molecule, space_group, z_value, z_prime = get_crystal_properties(structure_identifier)
    crystal_path = crystals_path + f'{structure_identifier}.pdb'
    print("Loading Crystal" + crystal_path)
    # gather structural information from the crystal file
    unit_cell = io.read(crystal_path)
    crystal_atoms = unit_cell.get_atomic_numbers()
    single_mol_atoms = crystal_atoms[:atoms_in_molecule]
    crystal_coordinates = unit_cell.positions
    cell_params = unit_cell.cell
    cell_lengths, cell_angles = unit_cell.cell.lengths(), unit_cell.cell.angles()
    T_fc = unit_cell.cell.array  # fractional to cartesian coordinate transform matrix

    # generate supercell
    supercell_coordinates = []
    for xs in range(cluster_size[0]):
        for ys in range(cluster_size[1]):
            for zs in range(cluster_size[2]):
                supercell_coordinates.extend(crystal_coordinates + T_fc[0] * xs + T_fc[1] * ys + T_fc[2] * zs)

    supercell_coordinates = np.asarray(supercell_coordinates)
    supercell_atoms = np.concatenate([crystal_atoms for _ in range(np.prod(cluster_size))])

    # adjust shape of the cluster
    if cluster_type == "supercell":
        cluster_size, supercell_atoms, supercell_coordinates =(
            build_supercell(T_fc, cell_lengths, cluster_size, crystal_atoms, crystal_coordinates,
                            min_lattice_length, supercell_atoms))
    elif cluster_type == "spherical":  # exclude molecules beyond some radial cutoff
        supercell_atoms, supercell_coordinates = (carve_spherical_cluster(
                atoms_in_molecule, cell_lengths, cluster_size, max_sphere_radius, single_mol_atoms,
                supercell_atoms, supercell_coordinates, z_value))

    if prep_crystal_in_melt:
        melt_inds, supercell_coordinates = (
            crystal_melt_reindexing(
                atoms_in_molecule, cluster_size, max_sphere_radius, supercell_coordinates, z_value))

    else:
        melt_inds = None

    if scramble_rate > 0:
        supercell_atoms, supercell_coordinates =\
            apply_scramble(
                atoms_in_molecule, scramble_rate, single_mol_atoms, supercell_atoms, supercell_coordinates)

    if gap_rate > 0:
        supercell_atoms, supercell_coordinates = (
            apply_gap(
                atoms_in_molecule, gap_rate, single_mol_atoms, supercell_atoms, supercell_coordinates))

    if defect_rate > 0:  # sub nicotinamides for benzamides
        if 'acridine' in crystal_path:
            assert False, "Cannot do defects for acridine yet."
        supercell_atoms, supercell_coordinates = (
            apply_defect(
                atoms_in_molecule, defect_rate, single_mol_atoms, supercell_atoms, supercell_coordinates))

    if periodic_structure:
        cell = T_fc * np.asarray(cluster_size)  # cell parameters are the same as the fractional->cartesian transition matrix (or sometimes its transpose)
    else:
        cell = (np.ptp(supercell_coordinates) + min_inter_cluster_distance) * np.eye(3) / 2

    supercell_coordinates += cell.sum(0) / 2 - supercell_coordinates.mean(0)

    cluster = Atoms(positions=supercell_coordinates, numbers=supercell_atoms, cell=cell)

    filename = f'{"_".join(structure_identifier.split("/"))}_{cluster_type}_{cluster_size}_defect={defect_rate}_vacancy={gap_rate}_disorder={scramble_rate}.xyz'
    io.write(filename, cluster)

    return filename, melt_inds


def crystal_melt_reindexing(atoms_in_molecule, cluster_size, max_sphere_radius, supercell_coordinates, z_value):
    # identify atoms in a sufficiently large sphere from the center
    # reindex the whole thing to put these in the first N rows
    num_mols = z_value * np.product(cluster_size)
    molwise_supercell_coordinates = supercell_coordinates.reshape(num_mols, atoms_in_molecule, 3)
    centroid = supercell_coordinates.mean(0)
    mol_centroids = molwise_supercell_coordinates.mean(1)
    dists = cdist(centroid[None, :], mol_centroids)[0, :]
    crystal_mol_inds = np.argwhere(dists < max_sphere_radius)[:, 0]
    melt_mol_inds = np.argwhere(dists > max_sphere_radius)[:, 0]
    # complete reindexing
    molwise_supercell_coordinates = np.concatenate([molwise_supercell_coordinates[crystal_mol_inds], molwise_supercell_coordinates[melt_mol_inds]], axis=0)
    # supercell atom indexing doesn't change
    supercell_coordinates = molwise_supercell_coordinates.reshape(int(len(molwise_supercell_coordinates) * atoms_in_molecule), 3)
    melt_inds = {'melt_start_ind': len(crystal_mol_inds) + 1,
                 'melt_end_ind': len(molwise_supercell_coordinates),
                 'crystal_start_ind': 1,
                 'crystal_end_ind': len(crystal_mol_inds)}
    melt_inds = dict2namespace(melt_inds)
    return melt_inds, supercell_coordinates


def build_supercell(T_fc, cell_lengths, cluster_size, crystal_atoms, crystal_coordinates, min_lattice_length, supercell_atoms):
    if min_lattice_length is not None:
        required_repeats = np.ceil(min_lattice_length / cell_lengths).astype(int)
        supercell_coordinates = []
        for xs in range(required_repeats[0]):
            for ys in range(required_repeats[1]):
                for zs in range(required_repeats[2]):
                    supercell_coordinates.extend(crystal_coordinates + T_fc[0] * xs + T_fc[1] * ys + T_fc[2] * zs)

        supercell_coordinates = np.asarray(supercell_coordinates)
        supercell_atoms = np.concatenate([crystal_atoms for _ in range(np.prod(required_repeats))])

        cluster_size = required_repeats

    else:
        pass
    return cluster_size, supercell_atoms, supercell_coordinates


def carve_spherical_cluster(atoms_in_molecule, cell_lengths, cluster_size, max_sphere_radius, single_mol_atoms, supercell_atoms, supercell_coordinates, z_value):
    num_mols = z_value * np.product(cluster_size)
    molwise_supercell_coordinates = supercell_coordinates.reshape(num_mols, atoms_in_molecule, 3)
    centroid = supercell_coordinates.mean(0)
    mol_centroids = molwise_supercell_coordinates.mean(1)
    dists = cdist(centroid[None, :], mol_centroids)[0, :]
    # find maximal spherical radius
    if max_sphere_radius is None:
        supercell_lengths = cell_lengths * cluster_size
        max_radius = min(supercell_lengths)
    else:
        max_radius = max_sphere_radius
    mols_to_keep = np.argwhere(dists < max_radius)[:, 0]
    keeper_molwise_coordinates = molwise_supercell_coordinates[mols_to_keep]
    # prune lonely molecules
    no_lonely = False
    loop_ind = 0
    while not no_lonely:
        loop_ind += 1
        # get centroid dists
        kept_molecule_centroids = keeper_molwise_coordinates.mean(1)
        kept_dists = cdist(kept_molecule_centroids, kept_molecule_centroids)

        # get coordination number
        coordination_number = np.zeros(len(kept_dists))
        max_mol_radius = np.amax(np.linalg.norm(molwise_supercell_coordinates[0] - molwise_supercell_coordinates[0].mean(0), axis=-1))
        inter_mol_dist = 2 * max_mol_radius  # 2 radii - empirically decided
        for ind in range(len(kept_dists)):
            coordination_number[ind] = np.sum(kept_dists[ind] < inter_mol_dist)

        # check if converged
        non_lonely_inds = np.argwhere(coordination_number > 2)[:, 0]  # 2 - also empirically decided
        if len(non_lonely_inds) == len(kept_molecule_centroids):
            no_lonely = True
        else:
            keeper_molwise_coordinates = keeper_molwise_coordinates[non_lonely_inds]  # update molwise coords
    # assign final coords & atomic numbers
    supercell_coordinates = keeper_molwise_coordinates.reshape(int(len(non_lonely_inds) * atoms_in_molecule), 3)
    supercell_atoms = np.concatenate([single_mol_atoms for _ in range(len(non_lonely_inds))])
    return supercell_atoms, supercell_coordinates


def apply_scramble(atoms_in_molecule, scramble_rate, single_mol_atoms, supercell_atoms, supercell_coordinates):
    num_mols = len(supercell_coordinates) // atoms_in_molecule
    num_defect_molecules = int(scramble_rate * num_mols)
    defect_molecule_indices = np.random.choice(np.arange(num_mols), size=num_defect_molecules)
    molwise_supercell_coordinates = supercell_coordinates.reshape(num_mols, atoms_in_molecule, 3)
    # apply defect
    defected_supercell_coordinates = []
    defected_supercell_atoms = []
    for i in range(len(molwise_supercell_coordinates)):
        if i in defect_molecule_indices:  # yes defect
            original_mol_coords = molwise_supercell_coordinates[i]
            random_rotation = Rotation.random().as_matrix()
            centroid = original_mol_coords.mean(0)
            rotated_mol_coords = np.inner(random_rotation, original_mol_coords - centroid).T + centroid
            defected_supercell_coordinates.append(rotated_mol_coords)
        else:  # no defect
            defected_supercell_coordinates.append(molwise_supercell_coordinates[i])
    supercell_atoms = np.concatenate([single_mol_atoms for _ in range(len(defected_supercell_coordinates))])
    supercell_coordinates = np.concatenate(defected_supercell_coordinates)
    return supercell_atoms, supercell_coordinates


def apply_gap(atoms_in_molecule, gap_rate, single_mol_atoms, supercell_atoms, supercell_coordinates):
    num_mols = len(supercell_coordinates) // atoms_in_molecule
    num_defect_molecules = int(gap_rate * num_mols)
    defect_molecule_indices = np.random.choice(np.arange(num_mols), size=num_defect_molecules)
    molwise_supercell_coordinates = supercell_coordinates.reshape(num_mols, atoms_in_molecule, 3)
    # apply defect
    defected_supercell_coordinates = []
    defected_supercell_atoms = []
    for i in range(len(molwise_supercell_coordinates)):
        if i in defect_molecule_indices:  # yes defect
            pass
        else:  # no defect
            defected_supercell_coordinates.append(molwise_supercell_coordinates[i])
    supercell_atoms = np.concatenate([single_mol_atoms for _ in range(len(defected_supercell_coordinates))])
    supercell_coordinates = np.concatenate(defected_supercell_coordinates)
    return supercell_atoms, supercell_coordinates


def apply_defect(atoms_in_molecule, defect_rate, single_mol_atoms, supercell_atoms, supercell_coordinates):
    benzamide_atoms = np.concatenate((single_mol_atoms, np.ones(1))).astype(int)  # add a hydrogen # TODO get this in a non-hardcoded fashion
    atom_switch_coord = 2
    benzamide_atoms[atom_switch_coord] = 6  # replace ring nitrogen by carbon
    BENZINE_C_H_BOND_LENGTH = 1.09  # angstrom
    # pick molecules to defect
    num_mols = len(supercell_coordinates) // atoms_in_molecule
    num_defect_molecules = int(defect_rate * num_mols)
    defect_molecule_indices = np.random.choice(np.arange(num_mols), size=num_defect_molecules)
    molwise_supercell_coordinates = supercell_coordinates.reshape(num_mols, atoms_in_molecule, 3)
    # apply defect
    defected_supercell_coordinates = []
    defected_supercell_atoms = []
    for i in range(len(molwise_supercell_coordinates)):
        if i in defect_molecule_indices:  # yes defect
            original_mol_coords = molwise_supercell_coordinates[i]

            # append a proton in the correct spot
            new_carbon_coord = original_mol_coords[atom_switch_coord]
            neighboring_carbon_inds = list(np.argwhere(cdist(new_carbon_coord[None, :], original_mol_coords)[0, :] < 1.45)[:, 0])
            neighboring_carbon_inds.remove(2)  # remove self
            neighbor_vectors = new_carbon_coord - original_mol_coords[neighboring_carbon_inds]

            # to project a trigonal planar proton, take the mean of the neighbors directions, and reverse it
            normed_neighbor_vectors = neighbor_vectors / np.linalg.norm(neighbor_vectors, axis=1)[:, None]
            proton_direction = normed_neighbor_vectors.mean(0)  # switch direction
            proton_vector = proton_direction / np.linalg.norm(proton_direction) * BENZINE_C_H_BOND_LENGTH
            proton_position = new_carbon_coord + proton_vector

            defect_mol_coordinates = np.concatenate((original_mol_coords, proton_position[None, :]))  # append proton position to end of list

            defected_supercell_coordinates.append(defect_mol_coordinates)
            defected_supercell_atoms.extend(benzamide_atoms)
        else:  # no defect
            defected_supercell_coordinates.append(molwise_supercell_coordinates[i])
            defected_supercell_atoms.extend(single_mol_atoms)
    supercell_coordinates = np.concatenate(defected_supercell_coordinates)
    supercell_atoms = np.asarray(defected_supercell_atoms)
    return supercell_atoms, supercell_coordinates

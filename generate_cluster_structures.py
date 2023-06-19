"""
script for converting unit cells to finite size crystal structures
options for cluster size and shape
options for introduction of defects
"""

import os
import sys

import numpy as np
from ase import Atoms
from ase.visualize import view
from ase import io
from scipy.spatial.distance import cdist, pdist
from scipy.spatial.transform import Rotation

# options
structure_identifier = "NICOAM07"
cluster_type = "supercell"  # "supercell" or "spherical"
max_sphere_radius = 20
cluster_size = [1, 1, 1]  # size of supercell to be generated [a,b,c]
defect_rate = 0  # fraction of molecules which will be switched from nicotinamide to benzamide
scramble_rate = 0  # fraction of molecules with scrambled orientations
gap_rate = 0  # fraction of molecules deleted
seed = 0

np.random.seed(seed=seed)

# move to working directory
os.chdir(r'C:\Users\mikem\crystals\clusters\cluster_structures\test_2')

# get original structure
if structure_identifier == "NICOAM07":
    atoms_in_molecule = 15
    space_group = "P21"
    z_value = 2
    crystal_path = r'C:\Users\mikem\crystals\clusters\Leslie\CrystalStructures\NICOAM07_renumbered_1x1x1.pdb'
elif structure_identifier == "NICOAM013":
    atoms_in_molecule = 15
    space_group = "P21/c"
    z_value = 2
    crystal_path = r'C:\Users\mikem\crystals\clusters\Leslie\CrystalStructures\NICOAM13_renumbered_1x1x1.pdb'
else:
    print("no such structure!")
    sys.exit()

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
    pass
elif cluster_type == "spherical":  # exclude molecules beyond some radial cutoff
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

    keeper_molecule_coordinates = molwise_supercell_coordinates[mols_to_keep].reshape(int(len(mols_to_keep) * atoms_in_molecule), 3)

    supercell_coordinates = keeper_molecule_coordinates
    supercell_atoms = np.concatenate([single_mol_atoms for _ in range(len(mols_to_keep))])

if scramble_rate > 0:
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

    supercell_coordinates = np.concatenate(defected_supercell_coordinates)

if gap_rate > 0:
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

    supercell_coordinates = np.concatenate(defected_supercell_coordinates)

if defect_rate > 0:  # sub nicotinamides for benzamides
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

cell = np.ptp(supercell_coordinates)*10 * np.eye(3)
supercell_coordinates += cell.sum(0) / 2
cluster = Atoms(positions=supercell_coordinates, numbers=supercell_atoms, cell=cell)  # cell = T_fc)

io.write('test_cluster.xyz', cluster)

aa = 1

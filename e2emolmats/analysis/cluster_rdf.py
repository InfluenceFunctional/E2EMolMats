import itertools
import torch

from e2emolmats.analysis.rdf_utils import parallel_compute_rdf_torch, asymmetric_radius_graph


def crystal_rdf(positions, symbols, mol_num_atoms, in_inds, out_inds, rrange=[0, 10], bins=100, elementwise=False, raw_density=False, atomwise=False, cpu_detach=True):
    """
    compute rdf for a single crystal
    """
    positions = torch.Tensor(positions)
    symbols = torch.tensor(symbols, dtype=torch.int64)
    in_inds = torch.tensor(in_inds, dtype=torch.int64)
    out_inds = torch.tensor(out_inds, dtype=torch.int64)
    batch = torch.zeros(len(positions), dtype=torch.int64)
    num_graphs = 1
    device = 'cpu'

    edges = asymmetric_radius_graph(
        positions,
        batch=batch,
        inside_inds=in_inds,
        convolve_inds=out_inds,
        r=max(rrange), max_num_neighbors=500, flow='source_to_target')

    # track which edges go with which crystals
    edge_in_crystal_number = batch[edges[0]]

    # compute all the dists
    dists = (positions[edges[0]] - positions[edges[1]]).pow(2).sum(dim=-1).sqrt()

    assert not (elementwise and atomwise)

    if elementwise:
        relevant_elements = [5, 6, 7, 8, 9, 15, 16, 17, 35]
        element_symbols = {5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 35: 'Br'}
        elements = [symbols[edges[0]], symbols[edges[1]]]
        rdfs_dict = {}
        rdfs_array = torch.zeros((num_graphs, int((len(relevant_elements) ** 2 + len(relevant_elements)) / 2), bins))
        ind = 0
        for i, element1 in enumerate(relevant_elements):
            for j, element2 in enumerate(relevant_elements):
                if j >= i:
                    rdfs_dict[ind] = element_symbols[element1] + ' to ' + element_symbols[element2]
                    rdfs_array[:, ind], rr = \
                        parallel_compute_rdf_torch(
                            [dists[(edge_in_crystal_number == n) * (elements[0] == element1) * (elements[1] == element2)]
                             for n in range(num_graphs)],
                            rrange=rrange, bins=bins, raw_density=raw_density)
                    ind += 1
        if cpu_detach:
            rdfs_array = rdfs_array.cpu().detach().numpy()
            rr = rr.cpu().detach().numpy()

        return rdfs_array, rr, rdfs_dict

    elif atomwise:  # generate atomwise indices which are shared between samples
        rdfs_array_list = []
        rdfs_dict_list = []
        all_atom_inds = []
        # abs_atom_inds = []
        for i in range(num_graphs):
            # # assume only that the order of atoms is patterned in all images
            # canonical_conformer_coords = positions[:mol_num_atoms]
            # centroid = canonical_conformer_coords.mean(0)
            # mol_dists = torch.linalg.norm(canonical_conformer_coords - centroid[None, :], dim=-1)
            # # todo if any dists are within 32 bit uncertainty, we have to break the symmetry somehow - can lead to errors in RDF comparisons
            # inds = torch.argsort(mol_dists)
            # all_atom_inds.append(inds.tile(int((batch == i).sum() // mol_num_atoms)))
            all_atom_inds.append(torch.arange(mol_num_atoms).tile(int((batch == i).sum() // mol_num_atoms)))  # assume our molecules are always in the same order

        atom_inds = torch.cat(all_atom_inds)

        atoms_in_edges = [atom_inds[edges[0]], atom_inds[edges[1]]]

        ind = 0
        for n in range(num_graphs):
            # this way is slightly faster than the above, and should scale well to larger batch sizes
            # all possible combinations of unique atoms on this graph
            atom_pairs = torch.Tensor(list(itertools.combinations(torch.arange(int(mol_num_atoms)), 2)))
            atom_pairs = torch.cat([atom_pairs, torch.arange(15).tile(2, 1).T], dim=0)  # above does not include self connections

            rdfs_dict_list.append(atom_pairs)  # record the pairs for reporting purposes

            in_crystal_inds = torch.where(edge_in_crystal_number == n)[0]  # atom indices in this crystal

            atom_locations = [[atoms_in_edges[0][in_crystal_inds] == m, atoms_in_edges[1][in_crystal_inds] == m] for m in range(int(atom_pairs.max()) + 1)]

            relevant_atoms_dists_list = [dists[
                                             in_crystal_inds[
                                                 torch.logical_and(
                                                     atom_locations[int(atom_pairs[m, 0])][0],
                                                     atom_locations[int(atom_pairs[m, 1])][1])
                                             ]
                                         ]
                                         for m in range(len(atom_pairs))]

            rdfs_array, rr = parallel_compute_rdf_torch(relevant_atoms_dists_list,
                                                        rrange=rrange, bins=bins, raw_density=raw_density)
            ind += 1
            if cpu_detach:
                rdfs_array = rdfs_array.cpu().detach().numpy()
                rr = rr.cpu().detach().numpy()

            rdfs_array_list.append(rdfs_array)

        return rdfs_array_list, rr, rdfs_dict_list
    else:  # all-to-all mode
        rdfs_array, rr = parallel_compute_rdf_torch([dists[edge_in_crystal_number == n] for n in range(num_graphs)], rrange=rrange, bins=bins, raw_density=raw_density)
        if cpu_detach:
            rdfs_array = rdfs_array.cpu().detach().numpy()
            rr = rr.cpu().detach().numpy()

        return rdfs_array, rr

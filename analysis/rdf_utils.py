from typing import Optional
import torch_cluster
import torch


def torch_ptp(tensor: torch.tensor):
    return torch.max(tensor) - torch.min(tensor)


def parallel_compute_rdf_torch(dists_list, raw_density=True, rrange=None, bins=None, remove_radial_scaling=False):
    """
    compute the radial distribution for a single fixed point
    dists: array of pairwise distances of nearby particles from the reference
    some batching for speed
    """
    hist_range = [0.5, 10] if rrange is None else rrange
    hist_bins = 100 if bins is None else bins

    if not raw_density:  # estimate the density from the distances
        rdf_density = torch.zeros(len(dists_list)).to(dists_list[0].device)
        for i in range(len(dists_list)):
            dists = dists_list[i]
            try:
                sorted_dists = torch.sort(dists)[0][:len(dists) // 2]  # we will use 1/2 the dist radius to avoid edge / cutoff effects
                volume = 4 / 3 * torch.pi * torch_ptp(sorted_dists) ** 3  # volume of a sphere #np.ptp(sorted_dists[:, 0]) * np.ptp(sorted_dists[:, 1]) * np.ptp(sorted_dists[:, 2])
                # number of particles divided by the volume
                rdf_density[i] = len(sorted_dists) / volume
            except:
                rdf_density[i] = 1
    else:
        rdf_density = torch.ones(len(dists_list), device=dists_list[0].device, dtype=torch.float32)

    hh_list = torch.stack([torch.histc(dists, min=hist_range[0], max=hist_range[1], bins=hist_bins) for dists in dists_list])
    rr = torch.linspace(hist_range[0], hist_range[1], hist_bins + 1).to(hh_list.device)
    if remove_radial_scaling:
        rdf = hh_list / rdf_density[:, None]  # un-smoothed radial density
    else:
        shell_volumes = (4 / 3) * torch.pi * ((rr[:-1] + torch.diff(rr)) ** 3 - rr[:-1] ** 3)  # volume of the shell at radius r+dr
        rdf = hh_list / shell_volumes[None, :] / rdf_density[:, None]  # un-smoothed radial density

    return rdf, (rr[:-1] + torch.diff(rr)).requires_grad_()  # rdf and x-axis


def radius(x: torch.Tensor, y: torch.Tensor, r: float,
           batch_x: Optional[torch.Tensor] = None,
           batch_y: Optional[torch.Tensor] = None,
           max_num_neighbors: int = 32,
           num_workers: int = 1) -> torch.Tensor:
    r"""Finds for each element in :obj:`y` all points in :obj:`x` within
    distance :obj:`r`.

    Args:
        x (Tensor): Node feature matrix
            :math:`\mathbf{X} \in \mathbb{R}^{N \times F}`.
        y (Tensor): Node feature matrix
            :math:`\mathbf{Y} \in \mathbb{R}^{M \times F}`.
        r (float): The radius.
        batch_x (LongTensor, optional): Batch vector
            :math:`\mathbf{b} \in {\{ 0, \ldots, B-1\}}^N`, which assigns each
            node to a specific example. :obj:`batch_x` needs to be sorted.
            (default: :obj:`None`)
        batch_y (LongTensor, optional): Batch vector
            :math:`\mathbf{b} \in {\{ 0, \ldots, B-1\}}^M`, which assigns each
            node to a specific example. :obj:`batch_y` needs to be sorted.
            (default: :obj:`None`)
        max_num_neighbors (int, optional): The maximum number of neighbors to
            return for each element in :obj:`y`.
            If the number of actual neighbors is greater than
            :obj:`max_num_neighbors`, returned neighbors are picked randomly.
            (default: :obj:`32`)
        num_workers (int): Number of workers to use for computation. Has no
            effect in case :obj:`batch_x` or :obj:`batch_y` is not
            :obj:`None`, or the input lies on the GPU. (default: :obj:`1`)

    .. code-block:: python

        import torch
        from torch_cluster import radius

        x = torch.Tensor([[-1, -1], [-1, 1], [1, -1], [1, 1]])
        batch_x = torch.tensor([0, 0, 0, 0])
        y = torch.Tensor([[-1, 0], [1, 0]])
        batch_y = torch.tensor([0, 0])
        assign_index = radius(x, y, 1.5, batch_x, batch_y)
    """

    x = x.view(-1, 1) if x.dim() == 1 else x
    y = y.view(-1, 1) if y.dim() == 1 else y
    x, y = x.contiguous(), y.contiguous()

    ptr_x: Optional[torch.Tensor] = None
    if batch_x is not None:
        assert x.size(0) == batch_x.numel()
        batch_size = int(batch_x.max()) + 1

        deg = x.new_zeros(batch_size, dtype=torch.long)
        deg.scatter_add_(0, batch_x, torch.ones_like(batch_x))

        ptr_x = deg.new_zeros(batch_size + 1)
        torch.cumsum(deg, 0, out=ptr_x[1:])

    ptr_y: Optional[torch.Tensor] = None
    if batch_y is not None:
        assert y.size(0) == batch_y.numel()
        batch_size = int(batch_y.max()) + 1

        deg = y.new_zeros(batch_size, dtype=torch.long)
        deg.scatter_add_(0, batch_y, torch.ones_like(batch_y))

        ptr_y = deg.new_zeros(batch_size + 1)
        torch.cumsum(deg, 0, out=ptr_y[1:])

    return torch.ops.torch_cluster.radius(x, y, ptr_x, ptr_y, r,
                                          max_num_neighbors, num_workers)


def asymmetric_radius_graph(x: torch.Tensor, r: float,
                            inside_inds: torch.Tensor, convolve_inds: torch.Tensor,
                            batch: torch.Tensor,
                            loop: bool = False,
                            max_num_neighbors: int = 32, flow: str = 'source_to_target',
                            num_workers: int = 1) -> torch.Tensor:
    r"""Computes graph edges to all points within a given distance.

    Args:
        x (Tensor): Node feature matrix
            :math:`\mathbf{X} \in \mathbb{R}^{N \times F}`.
        r (float): The radius.
        batch (LongTensor, optional): Batch vector
            :math:`\mathbf{b} \in {\{ 0, \ldots, B-1\}}^N`, which assigns each
            node to a specific example. :obj:`batch` needs to be sorted.
            (default: :obj:`None`)
        loop (bool, optional): If :obj:`True`, the graph will contain
            self-loops. (default: :obj:`False`)
        max_num_neighbors (int, optional): The maximum number of neighbors to
            return for each element.
            If the number of actual neighbors is greater than
            :obj:`max_num_neighbors`, returned neighbors are picked randomly.
            (default: :obj:`32`)
        flow (string, optional): The flow direction when used in combination
            with message passing (:obj:`"source_to_target"` or
            :obj:`"target_to_source"`). (default: :obj:`"source_to_target"`)
        num_workers (int): Number of workers to use for computation. Has no
            effect in case :obj:`batch` is not :obj:`None`, or the input lies
            on the GPU. (default: :obj:`1`)
        inside_inds (Tensor): original indices for the nodes in the y subgraph

    :rtype: :class:`LongTensor`

    .. code-block:: python

        import torch
        from torch_cluster import radius_graph

        x = torch.Tensor([[-1, -1], [-1, 1], [1, -1], [1, 1]])
        batch = torch.tensor([0, 0, 0, 0])
        edge_index = radius_graph(x, r=1.5, batch=batch, loop=False)
    """
    if convolve_inds is None:  # indexes of items within x to convolve against y
        convolve_inds = torch.arange(len(x))

    assert flow in ['source_to_target', 'target_to_source']
    if batch is not None:
        edge_index = radius(x[convolve_inds], x[inside_inds], r, batch[convolve_inds], batch[inside_inds],
                            max_num_neighbors if loop else max_num_neighbors + 1,
                            num_workers)
    else:
        edge_index = radius(x[convolve_inds], x[inside_inds], r, None, None,
                            max_num_neighbors if loop else max_num_neighbors + 1,
                            num_workers)

    target, source = edge_index[0], edge_index[1]

    # edge_index[1] = inside_inds[edge_index[1, :]] # reindex
    target = inside_inds[target]  # contains correct indexes
    source = convolve_inds[source]

    if flow == 'source_to_target':
        row, col = source, target
    else:
        row, col = target, source

    if not loop:  # now properly deletes self-loops
        mask = row != col
        row, col = row[mask], col[mask]

    return torch.stack([row, col], dim=0)

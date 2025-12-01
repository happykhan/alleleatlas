"""Build an MST over coarse supernodes derived from a k-NN graph and plot it.

Pipeline:
  1. Load and symmetrize the k-NN matrix.
  2. Within each connected component, apply single-linkage clustering
     (MST + cutoff ``coarse_threshold``) to define coarse supernodes.
  3. Pick a medoid for each supernode using shortest-path distances.
  4. Build a supernode graph over medoids using shortest-path distances and
     take its MST.
  5. Lay out the supernode MST via ``layout_mode``:
     - "topology": unweighted spring layout of the MST.
     - "mds": MDS on medoid-to-medoid distances, then a short unweighted
       spring pass to spread nodes.
  6. Draw nodes with log-scaled sizes (or a fixed size) and MST edges
     (optionally filtered by ``edge_display_threshold``) shrunk to meet the
     node circumferences.

Many legacy parameters are kept for API compatibility but are ignored in the
current simplified layout (e.g., supernode_layout/component_layout/global_use_weights).
"""
from __future__ import annotations

import numpy as np
from typing import Optional
import csv
import math
import scipy.sparse as sp
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from rich.console import Console

console = Console()

def load_knn(path: str) -> sp.csr_matrix:
    return sp.load_npz(path)

def single_linkage_components(sub: sp.spmatrix, threshold: Optional[float]):
    """Return connected components after single-linkage cut at ``threshold``.

    We build an MST for the subgraph, then union edges whose weight is
    <= threshold. This mirrors single-linkage clustering while working
    on sparse graphs without building a dense distance matrix.
    """
    if sub.shape[0] == 0:
        return []
    if threshold is None:
        return [list(range(sub.shape[0]))]
    try:
        Gsub = nx.from_scipy_sparse_array(sub, edge_attribute='weight')
    except Exception:
        Gsub = nx.from_scipy_sparse_matrix(sub, edge_attribute='weight')
    if Gsub.number_of_edges() == 0:
        return [[i] for i in range(sub.shape[0])]
    T = nx.minimum_spanning_tree(Gsub, weight='weight')
    from networkx.utils import UnionFind
    uf = UnionFind()
    for n in T.nodes():
        # access item to ensure it's added to the union-find structure
        _ = uf[n]
    for u, v, d in T.edges(data=True):
        if d.get('weight', 1.0) <= threshold:
            uf.union(u, v)
    comps = {}
    for n in T.nodes():
        root = uf[n]
        comps.setdefault(root, []).append(int(n))
    return list(comps.values())


def build_component_mst(M: sp.csr_matrix) -> nx.Graph:
    # Ensure symmetry and convert to sparse graph with finite weights
    A = (M + M.T) / 2
    # create networkx graph using edges with finite weights
    try:
        G = nx.from_scipy_sparse_array(A, edge_attribute='weight')
    except Exception:
        G = nx.from_scipy_sparse_matrix(A, edge_attribute='weight')
    return G


def layout_and_draw(G: nx.Graph, ax, node_size=30):
    # spring layout for component
    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos, node_size=node_size, ax=ax)
    nx.draw_networkx_edges(G, pos, ax=ax)


def plot_cluster_backbone(
    knn: str,
    out: str = 'cluster_backbone.png',
    coarse_threshold: float = 50,
    min_component_size: int = 10,
    max_supernodes: int = 800,
    max_local_nodes: int = 1000,
    component_summary_out: Optional[str] = None,
    supernode_size_base: float = 30.0,
    supernode_size_scale: float = 6.0,
    edge_display_threshold: Optional[float] = 600,
    edge_cap: Optional[float] = None,
    grid_scale: float = 1.2,
    color_edges_by_weight: bool = False,
    layout_mode: str = 'topology',
    edge_shrink_scale: float = 1.0,
    edge_shrink_pad: float = 1.0,
    draw_nodes: bool = True,
    fixed_node_size: Optional[float] = None,
) -> None:
    """
    Plot a supernode MST for a k-NN graph.

    Steps:
      1. Load and symmetrize the k-NN matrix.
      2. Inside each connected component, apply single-linkage clustering
         at ``coarse_threshold`` to define coarse supernodes.
      3. Pick a medoid per supernode via shortest-path sums.
      4. Build a medoid-to-medoid distance graph and take its MST.
      5. Lay out the MST with ``layout_mode``:
         - "topology": unweighted spring layout of the MST.
         - "mds": MDS on the medoid distance matrix, followed by a short
           unweighted spring pass.
      6. Draw nodes with log-scaled sizes (or fixed sizes) and MST edges
         (optionally filtered by ``edge_display_threshold``) trimmed to
         node circumferences.

    Parameters
    ----------
    knn : str
        Path to .npz k-NN sparse matrix.
    out : str
        Output image filename.
    coarse_threshold : float
        Single-linkage cutoff for forming supernodes.
    min_component_size : int
        Minimum supernode size to keep.
    max_supernodes : int
        Maximum number of supernodes to plot (keep largest).
    max_local_nodes : int
        Cap per-component size when computing medoids.
    component_summary_out : Optional[str]
        If set, write a TSV summary of supernodes.
    edge_display_threshold : Optional[float]
        Only draw MST edges with weight <= this (after capping).
    edge_cap : Optional[float]
        Cap supernode edge weights before layout/drawing.
    grid_scale : float
        Scale the figure size.
    layout_mode : {"topology", "mds"}
        Global layout strategy for the supernode MST.
    edge_shrink_scale, edge_shrink_pad : float
        Control how much edges are shrunk to meet node boundaries.
    draw_nodes : bool
        If False, draw only edges.
    fixed_node_size : Optional[float]
        If set, use this constant node size instead of log-scaling.

    Notes
    -----
    Legacy parameters (supernode_layout/component_layout/spring params/
    global_use_weights/local_use_weights/etc.) are retained for API
    compatibility but are ignored by the simplified layout.
    """
    layout_mode = (layout_mode or "topology").lower()
    if layout_mode not in {"topology", "mds"}:
        raise ValueError("layout_mode must be 'topology' or 'mds'")

    console.print(f"Loading k-NN graph from: {knn}")
    console.print(f"Collapse threshold: {coarse_threshold}" )
    console.print(f"Edge display threshold: {edge_display_threshold}")
    if edge_cap is not None:
        console.print(f"[blue]Capping supernode edge weights at {edge_cap}[/blue]")
    M = load_knn(knn)
    A = (M + M.T) / 2
    console.print("[blue]Using single-linkage (MST) to collapse nodes[/blue]")
    ncomp, labels = sp.csgraph.connected_components(A, directed=False)
    components = []
    for comp_id in range(ncomp):
        nodes = np.where(labels == comp_id)[0]
        if len(nodes) == 0:
            continue
        sub = A[nodes[:, None], nodes]
        clusters_local = single_linkage_components(sub, coarse_threshold)
        for comp in clusters_local:
            components.append([int(nodes[i]) for i in comp])

    components = [c for c in components if len(c) >= min_component_size]
    if len(components) > max_supernodes:
        components = sorted(components, key=len, reverse=True)[:max_supernodes]
        console.print(f"[yellow]Too many supernodes; keeping top {max_supernodes} by size[/yellow]")
    if not components:
        console.print("[yellow]No coarse components found for cluster backbone plot[/yellow]")
        return
    console.print(f"[blue]Supernodes after filtering:[/blue] {len(components)}")

    if component_summary_out:
        with open(component_summary_out, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(["supernode_id", "size", "coarse_threshold", "example_members"])
            for idx, comp in enumerate(components):
                members_preview = ','.join(str(x) for x in comp[:50])
                writer.writerow([idx, len(comp), coarse_threshold, members_preview])
        console.print(f"[green]Wrote component summary:[/green] {component_summary_out}")
    console.print(f"[blue]Supernodes after coarse collapse:[/blue] {len(components)}")

    # compute medoid per component using shortest-path distances within the component
    medoids = []
    comp_sizes = []
    for comp_nodes in components:
        comp_sizes.append(len(comp_nodes))
        # optionally cap local inset size for very large components
        if len(comp_nodes) > max_local_nodes:
            comp_nodes_local = comp_nodes[:max_local_nodes]
        else:
            comp_nodes_local = comp_nodes

        sub = A[comp_nodes_local, :][:, comp_nodes_local]

        # medoid by shortest-path sum
        dist_mat = sp.csgraph.dijkstra(sub, directed=False, unweighted=False)
        dist_sums = np.where(np.isfinite(dist_mat), dist_mat, 0).sum(axis=1)
        medoid_local_idx = int(np.argmin(dist_sums))
        medoids.append(comp_nodes_local[medoid_local_idx])

    # build supernode graph using medoid shortest paths
    medoid_indices = np.array(medoids, dtype=int)
    super_dists = sp.csgraph.dijkstra(A, indices=medoid_indices, directed=False, unweighted=False)[:, medoid_indices]
    super_dists = np.where(np.isfinite(super_dists), super_dists, super_dists.max(initial=1.0) * 2)
    super_edges_rows, super_edges_cols = np.triu_indices_from(super_dists, k=1)
    super_G = nx.Graph()
    super_G.add_nodes_from(range(len(components)))
    cap_count = 0
    for i, j in zip(super_edges_rows, super_edges_cols):
        w = float(super_dists[i, j])
        if edge_cap is not None and w > edge_cap:
            cap_count += 1
            w = edge_cap
        super_G.add_edge(int(i), int(j), weight=w)
    if edge_cap is not None:
        console.print(f"[blue]Capped {cap_count} supernode edges at {edge_cap}[/blue]")
    # draw edges from MST; always use weights to define topology, even in topology mode
    super_T = nx.minimum_spanning_tree(super_G, weight='weight')
    edges_to_draw = list(super_T.edges(data=True))
    if edge_display_threshold is not None:
        edges_to_draw = [(u, v, d) for u, v, d in edges_to_draw if d.get('weight', 0.0) <= edge_display_threshold]
        console.print(f"[blue]Drawing {len(edges_to_draw)} super-edges <= {edge_display_threshold}[/blue]")
    # --- simplified global layout: only 'topology' or 'mds' ---
    n_super = super_T.number_of_nodes()
    pos_super: dict[int, tuple[float, float]] = {}
    if n_super > 0:
        if layout_mode == "mds":
            from sklearn.manifold import MDS

            D = super_dists.copy()
            finite = np.isfinite(D)
            max_finite = D[finite].max() if finite.any() else 1.0
            D[~finite] = max_finite * 2.0

            mds = MDS(
                n_components=2,
                dissimilarity="precomputed",
                random_state=42,
            )
            coords = mds.fit_transform(D)
            pos_super = {i: (float(coords[i, 0]), float(coords[i, 1])) for i in range(n_super)}
            k = 1.0 / math.sqrt(n_super)
            pos_super = nx.spring_layout(
                super_T,
                pos=pos_super,
                seed=42,
                weight=None,
                iterations=100,
                k=k,
            )
        else:
            k = 1.5 / math.sqrt(n_super) if n_super > 0 else None
            pos_super = nx.spring_layout(
                super_T,
                seed=42,
                weight=None,
                iterations=300,
                k=k,
            )
    # ensure we have positions for all nodes and they are finite; fallback to spring if needed
    def _has_nan(ps):
        arr = np.array(list(ps.values())) if ps else np.array([])
        return arr.size > 0 and (np.isnan(arr).any() or np.isinf(arr).any())

    if len(pos_super) < super_T.number_of_nodes() or _has_nan(pos_super):
        console.print("[yellow]Position fallback: using spring layout for missing/empty positions[/yellow]")
        pos_super = nx.spring_layout(super_T, seed=42, weight=None)
    if _has_nan(pos_super):
        console.print("[yellow]Position fallback: forcing simple grid positions[/yellow]")
        pos_super = {n: (float(i), float(j)) for i, n in enumerate(super_T.nodes()) for j in [0]}
    # debug bounds
    if pos_super:
        coords = np.array(list(pos_super.values()))
        try:
            console.print(
                f"Supernode pos range x=({coords[:,0].min():.2f},{coords[:,0].max():.2f}) "
                f"y=({coords[:,1].min():.2f},{coords[:,1].max():.2f}) nodes={len(pos_super)} edges_shown={len(edges_to_draw)}",
                markup=False,
            )
        except Exception:
            pass

    # plot
    fig, ax = plt.subplots(figsize=(10 * grid_scale, 8 * grid_scale))
    colorbar_handle = None
    if fixed_node_size is not None:
        node_sizes = [fixed_node_size for _ in comp_sizes]
    else:
        sizes = np.array(comp_sizes, dtype=float)
        log_sizes = np.log1p(sizes)
        log_ptp = np.ptp(log_sizes)
        if log_ptp > 0:
            norm = (log_sizes - log_sizes.min()) / log_ptp
        else:
            norm = np.zeros_like(log_sizes)
        node_sizes = supernode_size_base + supernode_size_scale * (1.0 + 4.0 * norm)
        node_sizes = node_sizes.tolist()
    # draw supernode MST (filtered edges) with distance-colored edges
    if edges_to_draw:
        # per-node radii in points for shrinking edges to node circumference (+pad)
        radii_pts = {n: edge_shrink_scale * math.sqrt(sz / math.pi) + edge_shrink_pad for n, sz in zip(super_T.nodes(), node_sizes)}
        weights_draw = [d.get('weight', 1.0) for _, _, d in edges_to_draw]
        # establish data-unit radii based on axis extents and figure size
        xs = np.array([p[0] for p in pos_super.values()])
        ys = np.array([p[1] for p in pos_super.values()])
        xr = xs.max() - xs.min() if xs.size else 1.0
        yr = ys.max() - ys.min() if ys.size else 1.0
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        px_w = bbox.width * fig.dpi
        px_h = bbox.height * fig.dpi
        data_per_px_x = xr / px_w if px_w > 0 else 1.0
        data_per_px_y = yr / px_h if px_h > 0 else 1.0
        radii_data = {n: (r * data_per_px_x, r * data_per_px_y) for n, r in radii_pts.items()}

        if edge_display_threshold is not None and edge_display_threshold > 0:
            wmin_norm = coarse_threshold if coarse_threshold is not None else 0.0
            wmax_norm = edge_display_threshold
        else:
            wmin_norm = coarse_threshold if coarse_threshold is not None else 0.0
            wmax = max(weights_draw) if len(weights_draw) > 0 else 1.0
            wmax_norm = wmax if wmax > 0 else 1.0
        norm = plt.Normalize(vmin=wmin_norm, vmax=wmax_norm)
        cmap = mcolors.LinearSegmentedColormap.from_list("orange_seafoam", ["#f97306", "#38b2ac"])
        for (u, v, d) in edges_to_draw:
            p1 = np.array(pos_super[int(u)])
            p2 = np.array(pos_super[int(v)])
            vec = p2 - p1
            dist = np.linalg.norm(vec)
            w = float(d.get('weight', 1.0))
            color = cmap(norm(w)) if color_edges_by_weight else 'black'
            if dist > 1e-9:
                unit = vec / dist
                rx_u, ry_u = radii_data.get(int(u), (0.0, 0.0))
                rx_v, ry_v = radii_data.get(int(v), (0.0, 0.0))
                p1_off = p1 + unit * rx_u
                p2_off = p2 - unit * rx_v
            else:
                p1_off = p1
                p2_off = p2
            arrow = plt.matplotlib.patches.FancyArrowPatch(
                p1_off,
                p2_off,
                connectionstyle="arc3",
                linewidth=1.2,
                alpha=0.7,
                color=color,
                arrowstyle='-',
            )
            ax.add_patch(arrow)
        if color_edges_by_weight:
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cb_ax = ax.inset_axes([1.03, 0.2, 0.025, 0.6])
            colorbar_handle = fig.colorbar(sm, cax=cb_ax)
            colorbar_handle.set_label(f"Edge weight (clamped {wmin_norm:.2f}-{wmax_norm:.2f})")
    else:
        console.print("[yellow]No super-edges to draw after filtering[/yellow]")
    # if positions ended up empty, fall back to a simple spring layout
    if not pos_super:
        console.print("[yellow]Fallback: using spring layout for supernodes[/yellow]")
        pos_super = nx.spring_layout(super_T, seed=42)
    if draw_nodes:
        nx.draw_networkx_nodes(super_T, pos_super, node_size=node_sizes, node_color='tab:blue', alpha=0.9, ax=ax)
    # node size legend with nicer rounded labels and offset to the side
    if draw_nodes:
        try:
            def nice_round(val: float) -> int:
                if val <= 0:
                    return 1
                import math
                mag = 10 ** math.floor(math.log10(val))
                mant = val / mag
                for m in [1, 2, 5, 10]:
                    if mant <= m:
                        return int(m * mag)
                return int(10 * mag)

            raw_sizes = np.array(comp_sizes, dtype=float)
            if raw_sizes.size > 0:
                candidates = [raw_sizes.min(), np.quantile(raw_sizes, 0.5), raw_sizes.max()]
                rounded = sorted({nice_round(c) for c in candidates if c > 0}, reverse=True)[:3]
                handles = []
                log_sizes_all = np.log1p(raw_sizes)
                log_min = log_sizes_all.min()
                log_max = log_sizes_all.max()
                log_ptp = log_max - log_min if log_max > log_min else 1.0
                for s in rounded:
                    ls = np.log1p(float(s))
                    norm = (ls - log_min) / log_ptp if log_ptp > 0 else 0.0
                    bubble_size = supernode_size_base + supernode_size_scale * (1.0 + 4.0 * norm)
                    handles.append(plt.scatter([], [], s=bubble_size, color='tab:blue', alpha=0.9, label=f"{s} nodes"))
                if handles:
                    ax.legend(
                        handles=handles,
                        title="Supernode size",
                        loc="lower left",
                        bbox_to_anchor=(1.08, 0.02),
                        scatterpoints=1,
                        borderpad=1.3,
                        frameon=True,
                        labelspacing=1.0,
                    )
        except Exception:
            pass
    # autoscale to data and add padding to keep nodes in view even with legends
    try:
        xs = np.array([p[0] for p in pos_super.values()])
        ys = np.array([p[1] for p in pos_super.values()])
        if xs.size > 0 and ys.size > 0:
            pad_x = (xs.max() - xs.min()) * 0.1 + 1e-3
            pad_y = (ys.max() - ys.min()) * 0.1 + 1e-3
            ax.set_xlim(xs.min() - pad_x, xs.max() + pad_x)
            ax.set_ylim(ys.min() - pad_y, ys.max() + pad_y)
    except Exception:
        pass


    ax.axis('off')
    plt.tight_layout()
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.close(fig)
    console.print(f"[green]Wrote cluster backbone plot:[/green] {out}")

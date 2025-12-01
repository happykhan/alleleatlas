"""Hull-based supernode visualization for k-NN graphs."""

from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.sparse as sp
from matplotlib.patches import Circle, Polygon
from rich.console import Console
from scipy.spatial import ConvexHull

from alleleatlas.mst_and_plot import load_knn, single_linkage_components

console = Console()


def _symmetrize(M: sp.spmatrix) -> sp.spmatrix:
    return (M + M.T) / 2


def _build_supernodes(
    A: sp.spmatrix,
    coarse_threshold: float,
    min_component_size: int,
    max_supernodes: int,
) -> List[List[int]]:
    ncomp, labels = sp.csgraph.connected_components(A, directed=False)
    components: List[List[int]] = []
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
    return components


def _medoids_and_mst(
    A: sp.spmatrix,
    components: List[List[int]],
    max_local_nodes: int,
    edge_cap: Optional[float],
) -> Tuple[np.ndarray, List[int], nx.Graph, float, np.ndarray]:
    if not components:
        return np.array([], dtype=int), [], nx.Graph(), 1.0, np.zeros((0, 0), dtype=float)

    medoids: List[int] = []
    comp_sizes: List[int] = []
    for comp_nodes in components:
        comp_sizes.append(len(comp_nodes))
        comp_nodes_local = comp_nodes[:max_local_nodes] if len(comp_nodes) > max_local_nodes else comp_nodes
        sub = A[comp_nodes_local, :][:, comp_nodes_local]
        dist_mat = sp.csgraph.dijkstra(sub, directed=False, unweighted=False)
        dist_sums = np.where(np.isfinite(dist_mat), dist_mat, 0).sum(axis=1)
        medoid_local_idx = int(np.argmin(dist_sums))
        medoids.append(comp_nodes_local[medoid_local_idx])

    medoid_indices = np.array(medoids, dtype=int)
    super_dists = sp.csgraph.dijkstra(A, indices=medoid_indices, directed=False, unweighted=False)[:, medoid_indices]
    super_dists = np.asarray(super_dists, dtype=float)
    finite_mask = np.isfinite(super_dists)
    if finite_mask.any():
        max_finite = float(super_dists[finite_mask].max())
    else:
        max_finite = 1.0
    fill_value = max(max_finite * 2.0, 1.0)
    super_dists = np.where(finite_mask, super_dists, fill_value)

    n_super = len(components)
    super_G = nx.Graph()
    super_G.add_nodes_from(range(n_super))

    rows, cols = np.triu_indices(n_super, k=1)
    cap_count = 0
    for i, j in zip(rows, cols):
        w = float(super_dists[i, j])
        if edge_cap is not None and w > edge_cap:
            w = edge_cap
            cap_count += 1
        super_G.add_edge(int(i), int(j), weight=w)
    if edge_cap is not None and cap_count > 0:
        console.print(f"[blue]Capped {cap_count} supernode edges at {edge_cap}[/blue]")

    super_T = nx.minimum_spanning_tree(super_G, weight="weight")
    return medoid_indices, comp_sizes, super_T, max_finite, super_dists


def _layout_positions(
    super_T: nx.Graph,
    super_dists: np.ndarray,
    layout_mode: str,
    random_seed: int,
) -> Dict[int, Tuple[float, float]]:
    n_super = super_T.number_of_nodes()
    if n_super == 0:
        return {}

    layout_mode = (layout_mode or "mds").lower()
    pos_super: Dict[int, Tuple[float, float]] = {}
    if layout_mode not in {"mds", "topology"}:
        layout_mode = "mds"

    if layout_mode == "mds":
        try:
            from sklearn.manifold import MDS

            D = np.asarray(super_dists, dtype=float)
            mds = MDS(
                n_components=2,
                dissimilarity="precomputed",
                random_state=random_seed,
            )
            coords = mds.fit_transform(D)
            pos_super = {i: (float(coords[i, 0]), float(coords[i, 1])) for i in range(n_super)}

            k = 1.0 / math.sqrt(n_super) if n_super > 0 else None
            pos_super = nx.spring_layout(
                super_T,
                pos=pos_super,
                seed=random_seed,
                weight=None,
                iterations=20,
                k=k,
            )
        except Exception as e:
            console.print(f"[yellow]MDS failed ({e}); falling back to spring layout[/yellow]")
            layout_mode = "topology"

    if layout_mode == "topology":
        k = 1.5 / math.sqrt(n_super) if n_super > 0 else None
        pos_super = nx.spring_layout(super_T, seed=random_seed, weight=None, iterations=2000, k=k)

    def _has_nan(ps: Dict[int, Tuple[float, float]]) -> bool:
        if not ps:
            return False
        arr = np.array(list(ps.values()))
        return np.isnan(arr).any() or np.isinf(arr).any()

    if len(pos_super) < n_super or _has_nan(pos_super):
        console.print("[yellow]Position fallback: using spring layout[/yellow]")
        pos_super = nx.spring_layout(super_T, seed=random_seed, weight=None)

    if _has_nan(pos_super):
        console.print("[yellow]Position fallback: forcing simple grid positions[/yellow]")
        pos_super = {n: (float(i), 0.0) for i, n in enumerate(super_T.nodes())}

    coords_arr = np.array(list(pos_super.values()))
    if coords_arr.size > 0:
        center = coords_arr.mean(axis=0)
        coords_arr = coords_arr - center
        max_span = max(np.ptp(coords_arr, axis=0).max(), 1e-6)
        coords_arr = coords_arr / max_span
        for idx, node in enumerate(pos_super.keys()):
            pos_super[node] = (float(coords_arr[idx, 0]), float(coords_arr[idx, 1]))

    return pos_super


def _draw_hulls(
    ax: plt.Axes,
    pos_super: Dict[int, Tuple[float, float]],
    super_T: nx.Graph,
    hull_threshold: float,
    min_nodes_in_hull: int = 2,
    cmap_name: str = "tab20",
    face_alpha: float = 0.15,
) -> None:
    H = nx.Graph()
    H.add_nodes_from(super_T.nodes())
    for u, v, d in super_T.edges(data=True):
        w = d.get("weight", 1.0)
        if np.isfinite(w) and w <= hull_threshold:
            H.add_edge(u, v, weight=w)

    comps = [list(c) for c in nx.connected_components(H)]
    comps = [c for c in comps if len(c) >= min_nodes_in_hull]
    if not comps:
        console.print("[yellow]No hull components at this threshold[/yellow]")
        return

    cmap = cm.get_cmap(cmap_name, len(comps))
    for idx, comp in enumerate(comps):
        pts = np.array([pos_super[n] for n in comp if n in pos_super])
        if pts.shape[0] == 0:
            continue
        color = cmap(idx)

        if pts.shape[0] == 1:
            x, y = pts[0]
            circ = Circle((x, y), radius=0.05, facecolor=color, alpha=face_alpha, edgecolor=color, linewidth=1.0)
            ax.add_patch(circ)
        elif pts.shape[0] == 2:
            p1, p2 = pts
            vec = p2 - p1
            norm = np.linalg.norm(vec)
            if norm == 0:
                x, y = p1
                circ = Circle((x, y), radius=0.05, facecolor=color, alpha=face_alpha, edgecolor=color, linewidth=1.0)
                ax.add_patch(circ)
            else:
                unit = vec / norm
                ortho = np.array([-unit[1], unit[0]])
                pad = 0.03
                poly_pts = np.vstack(
                    [
                        p1 - pad * ortho,
                        p2 - pad * ortho,
                        p2 + pad * ortho,
                        p1 + pad * ortho,
                    ]
                )
                poly = Polygon(poly_pts, closed=True, facecolor=color, alpha=face_alpha, edgecolor=color, linewidth=1.0)
                ax.add_patch(poly)
        else:
            hull = ConvexHull(pts)
            hull_pts = pts[hull.vertices]
            poly = Polygon(hull_pts, closed=True, facecolor=color, alpha=face_alpha, edgecolor=color, linewidth=1.0)
            ax.add_patch(poly)

        cx, cy = pts.mean(axis=0)
        ax.text(cx, cy, str(len(comp)), ha="center", va="center", fontsize=6, color="black", alpha=0.9)


def plot_knn_hulls(
    knn: str,
    out: str = "cluster_hulls.png",
    coarse_threshold: float = 50.0,
    min_component_size: int = 10,
    max_supernodes: int = 800,
    max_local_nodes: int = 1000,
    edge_cap: Optional[float] = None,
    edge_display_threshold: Optional[float] = 600.0,
    hull_threshold: Optional[float] = None,
    min_hull_nodes: int = 5,
    layout_mode: str = "mds",
    grid_scale: float = 1.2,
    draw_edges: bool = False,
    draw_nodes: bool = True,
    fixed_node_size: Optional[float] = None,
    random_seed: int = 42,
) -> None:
    """Plot supernodes with hulls grouping short MST edges."""
    layout_mode = (layout_mode or "mds").lower()
    if layout_mode not in {"mds", "topology"}:
        raise ValueError("layout_mode must be 'mds' or 'topology'")

    console.print(f"[blue]Loading k-NN graph for hull plot from:[/blue] {knn}")
    console.print(f"Coarse collapse threshold: {coarse_threshold}")
    console.print(f"Edge display threshold: {edge_display_threshold}")

    M = load_knn(knn)
    A = _symmetrize(M)

    components = _build_supernodes(
        A,
        coarse_threshold=coarse_threshold,
        min_component_size=min_component_size,
        max_supernodes=max_supernodes,
    )
    if not components:
        console.print("[yellow]No coarse components found; nothing to plot[/yellow]")
        return
    console.print(f"[blue]Supernodes after filtering:[/blue] {len(components)}")

    medoids, comp_sizes, super_T, max_finite, super_dists = _medoids_and_mst(
        A,
        components=components,
        max_local_nodes=max_local_nodes,
        edge_cap=edge_cap,
    )

    edges_to_draw = list(super_T.edges(data=True))
    if edge_display_threshold is not None:
        kept_edges = [(u, v, d) for u, v, d in edges_to_draw if d.get("weight", 0.0) <= edge_display_threshold]
    else:
        kept_edges = edges_to_draw
    console.print(
        f"[blue]MST edges:[/blue] {len(edges_to_draw)}, "
        f"edges <= threshold for drawing/hulls: {len(kept_edges)}"
    )

    pos_super = _layout_positions(super_T, super_dists, layout_mode=layout_mode, random_seed=random_seed)

    fig, ax = plt.subplots(figsize=(10 * grid_scale, 8 * grid_scale))

    if hull_threshold is None:
        if edge_display_threshold is not None:
            hull_threshold = edge_display_threshold * 1.2
        else:
            hull_threshold = max(max_finite, 1.0)
    console.print(f"Using hull_threshold = {hull_threshold}")
    _draw_hulls(
        ax,
        pos_super=pos_super,
        super_T=super_T,
        hull_threshold=hull_threshold,
        min_nodes_in_hull=max(2, int(min_hull_nodes)),
    )

    if draw_edges and kept_edges:
        for (u, v, d) in kept_edges:
            p1 = np.array(pos_super[int(u)])
            p2 = np.array(pos_super[int(v)])
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], linewidth=0.5, color="black", alpha=0.25, zorder=1)

    if draw_nodes:
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
            base = 40.0
            scale = 200.0
            node_sizes = (base + scale * norm).tolist()

        nx.draw_networkx_nodes(
            super_T,
            pos_super,
            node_size=node_sizes,
            node_color="tab:blue",
            alpha=0.9,
            ax=ax,
            linewidths=0.5,
            edgecolors="black",
        )

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

    ax.axis("off")
    plt.tight_layout()
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    console.print(f"[green]Wrote hull plot:[/green] {out}")

"""Diagnostic plots for k-NN clustering and supernode MSTs."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.sparse as sp
from networkx.utils import UnionFind
from rich.console import Console

from alleleatlas.mst_and_plot import load_knn, single_linkage_components

console = Console()


def _save_fig(fig, path: Path) -> str:
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return str(path)


def _symmetrize(M: sp.spmatrix) -> sp.spmatrix:
    return (M + M.T) / 2


def _supernode_components(
    A: sp.spmatrix,
    coarse_threshold: float,
    min_component_size: int,
    max_supernodes: Optional[int],
) -> Tuple[List[List[int]], np.ndarray]:
    ncomp, labels = sp.csgraph.connected_components(A, directed=False)
    components: List[List[int]] = []
    for comp_id in range(ncomp):
        nodes = np.where(labels == comp_id)[0]
        if nodes.size == 0:
            continue
        sub = A[nodes[:, None], nodes]
        clusters_local = single_linkage_components(sub, coarse_threshold)
        for comp in clusters_local:
            components.append([int(nodes[i]) for i in comp])

    components = [c for c in components if len(c) >= min_component_size]
    if max_supernodes is not None and len(components) > max_supernodes:
        components = sorted(components, key=len, reverse=True)[:max_supernodes]
        console.print(f"[yellow]Too many supernodes; keeping top {max_supernodes} by size for diagnostics[/yellow]")

    comp_sizes = np.array([len(c) for c in components], dtype=int)
    return components, comp_sizes


def _medoids_and_super_mst(
    A: sp.spmatrix,
    components: List[List[int]],
    max_local_nodes: int,
    edge_cap: Optional[float],
) -> Tuple[np.ndarray, nx.Graph]:
    if not components:
        return np.array([], dtype=int), nx.Graph()

    medoids: List[int] = []
    for comp_nodes in components:
        if len(comp_nodes) > max_local_nodes:
            comp_nodes_local = comp_nodes[:max_local_nodes]
        else:
            comp_nodes_local = comp_nodes

        sub = A[comp_nodes_local, :][:, comp_nodes_local]
        dist_mat = sp.csgraph.dijkstra(sub, directed=False, unweighted=False)
        dist_sums = np.where(np.isfinite(dist_mat), dist_mat, 0).sum(axis=1)
        medoid_local_idx = int(np.argmin(dist_sums))
        medoids.append(comp_nodes_local[medoid_local_idx])

    medoid_indices = np.array(medoids, dtype=int)
    super_dists = sp.csgraph.dijkstra(A, indices=medoid_indices, directed=False, unweighted=False)[:, medoid_indices]
    finite = np.isfinite(super_dists)
    if finite.any():
        max_finite = super_dists[finite].max()
        fill_value = max(max_finite * 2, 1.0)
    else:
        fill_value = 1.0
    super_dists = np.where(finite, super_dists, fill_value)

    rows, cols = np.triu_indices_from(super_dists, k=1)
    super_G = nx.Graph()
    super_G.add_nodes_from(range(len(components)))
    for i, j in zip(rows, cols):
        w = float(super_dists[i, j])
        if edge_cap is not None and w > edge_cap:
            w = edge_cap
        super_G.add_edge(int(i), int(j), weight=w)

    super_T = nx.minimum_spanning_tree(super_G, weight="weight")
    return medoid_indices, super_T


def _mst_edges_sorted(A: sp.spmatrix) -> np.ndarray:
    T = sp.csgraph.minimum_spanning_tree(A)
    coo = T.tocoo()
    edges = np.vstack([coo.data, coo.row, coo.col]).T
    if edges.size:
        edges = edges[np.isfinite(edges[:, 0])]
        edges = edges[edges[:, 0].argsort()] if edges.size else edges
    else:
        edges = edges.reshape(0, 3)
    return edges


def _counts_vs_threshold(
    edges_sorted: np.ndarray,
    thresholds: Sequence[float],
    n: int,
    min_component_size: int,
) -> List[int]:
    uf = UnionFind(range(n))
    counts: List[int] = []
    edge_idx = 0
    for thr in thresholds:
        if not np.isfinite(thr):
            counts.append(0)
            continue
        while edge_idx < len(edges_sorted) and edges_sorted[edge_idx, 0] <= thr:
            _, u, v = edges_sorted[edge_idx]
            uf.union(int(u), int(v))
            edge_idx += 1
        root_to_size: Dict[int, int] = {}
        for i in range(n):
            r = uf[i]
            root_to_size[r] = root_to_size.get(r, 0) + 1
        counts.append(sum(1 for s in root_to_size.values() if s >= min_component_size))
    return counts


def run_cluster_diagnostics(
    knn: str,
    *,
    coarse_threshold: float = 50,
    min_component_size: int = 10,
    max_supernodes: Optional[int] = None,
    max_local_nodes: int = 1000,
    edge_display_threshold: Optional[float] = 600,
    thresholds_sweep: Optional[Iterable[float]] = None,
    sweep_points: int = 25,
    out_prefix: str = "cluster_diag",
    edge_cap: Optional[float] = None,
) -> Dict[str, str]:
    """Generate diagnostic plots to tune clustering thresholds and layouts."""
    out_prefix_path = Path(out_prefix)
    out_prefix_path.parent.mkdir(parents=True, exist_ok=True)

    console.print(f"[blue]Running cluster diagnostics on:[/blue] {knn}")
    M = load_knn(knn)
    A = _symmetrize(M)
    n = A.shape[0]

    components, comp_sizes = _supernode_components(
        A,
        coarse_threshold=coarse_threshold,
        min_component_size=min_component_size,
        max_supernodes=max_supernodes,
    )
    paths: Dict[str, str] = {}
    if comp_sizes.size == 0:
        console.print("[yellow]No supernodes found; skipping diagnostics[/yellow]")
        return paths

    # 1) Histogram of supernode sizes
    bins_sizes = min(75, max(10, int(np.sqrt(comp_sizes.size))))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(comp_sizes, bins=bins_sizes, color="#3182bd", alpha=0.85)
    ax.set_xlabel("Supernode size (#profiles)")
    ax.set_ylabel("Count")
    if comp_sizes.max(initial=0) > 20:
        ax.set_yscale("log")
    ax.set_title("Distribution of supernode sizes")
    paths["supernode_size_hist"] = _save_fig(fig, out_prefix_path.with_name(out_prefix_path.name + "_supernode_size_hist.png"))

    # Build supernode MST once for downstream plots
    medoids, super_T = _medoids_and_super_mst(
        A,
        components=components,
        max_local_nodes=max_local_nodes,
        edge_cap=edge_cap,
    )
    mst_weights = np.array([d.get("weight", 1.0) for _, _, d in super_T.edges(data=True)], dtype=float)
    mst_weights = mst_weights[np.isfinite(mst_weights)]

    # 2) Histogram of MST edge weights
    if mst_weights.size > 0:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(mst_weights, bins=50, color="#31a354", alpha=0.85)
        ax.set_xlabel("Supernode MST edge weight")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of supernode MST edge weights")
        paths["mst_weight_hist"] = _save_fig(fig, out_prefix_path.with_name(out_prefix_path.name + "_supernode_mst_weights.png"))
    else:
        console.print("[yellow]Supernode MST has no edges; skipping MST weight histogram[/yellow]")

    # 3) Kept vs discarded edges relative to edge_display_threshold
    if edge_display_threshold is not None and mst_weights.size > 0:
        kept = mst_weights[mst_weights <= edge_display_threshold]
        discarded = mst_weights[mst_weights > edge_display_threshold]
        fig, ax = plt.subplots(figsize=(6, 4))
        bins = np.histogram_bin_edges(mst_weights, bins=50)
        if kept.size:
            ax.hist(kept, bins=bins, alpha=0.65, label="kept", color="#3182bd")
        if discarded.size:
            ax.hist(discarded, bins=bins, alpha=0.65, label="discarded", color="#de2d26")
        ax.axvline(edge_display_threshold, color="red", linestyle="--", label="threshold")
        ax.set_xlabel("Supernode MST edge weight")
        ax.set_ylabel("Count")
        ax.set_title("Kept vs discarded supernode MST edges")
        ax.legend()
        paths["kept_vs_discarded"] = _save_fig(fig, out_prefix_path.with_name(out_prefix_path.name + "_mst_edges_kept_vs_discarded.png"))

    # 4) Scatter: supernode size vs mean incident MST edge weight
    if mst_weights.size > 0 and len(super_T.nodes) == len(comp_sizes):
        mean_w = []
        sizes = []
        for node in super_T.nodes():
            weights = [d.get("weight", 1.0) for _, _, d in super_T.edges(node, data=True)]
            if weights:
                w_finite = [float(w) for w in weights if np.isfinite(w)]
                if not w_finite:
                    continue
                mean_w.append(float(np.mean(w_finite)))
                sizes.append(int(comp_sizes[int(node)]))
        if mean_w:
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.scatter(sizes, mean_w, s=12, alpha=0.55, color="#756bb1")
            ax.set_xscale("log")
            ax.set_xlabel("Supernode size (#profiles, log)")
            ax.set_ylabel("Mean MST edge weight to neighbours")
            ax.set_title("Supernode size vs mean MST neighbour distance")
            paths["size_vs_mean_edge"] = _save_fig(fig, out_prefix_path.with_name(out_prefix_path.name + "_size_vs_mean_mst_weight.png"))

    # 5) Sweep coarse_threshold: #supernodes vs threshold
    edges_sorted = _mst_edges_sorted(A)
    if thresholds_sweep is not None:
        thresholds_list = sorted({float(t) for t in thresholds_sweep if np.isfinite(t)})
    elif edges_sorted.size:
        weights_all = edges_sorted[:, 0]
        weights_all = weights_all[np.isfinite(weights_all)]
        if weights_all.size:
            qs = np.linspace(0.0, 1.0, sweep_points)
            thresholds_list = list(np.unique(np.quantile(weights_all, qs)))
        else:
            thresholds_list = [coarse_threshold]
    else:
        thresholds_list = [coarse_threshold]

    counts = _counts_vs_threshold(
        edges_sorted,
        thresholds_list,
        n=n,
        min_component_size=min_component_size,
    )
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(thresholds_list, counts, marker="o", markersize=4, linewidth=1.2)
    ax.set_xlabel("Single-linkage threshold")
    ax.set_ylabel("# supernodes (>= min_component_size)")
    ax.set_title("Supernodes vs single-linkage threshold")
    ax.grid(True, alpha=0.3)
    paths["supernodes_vs_threshold"] = _save_fig(fig, out_prefix_path.with_name(out_prefix_path.name + "_supernodes_vs_threshold.png"))

    console.print(f"[green]Cluster diagnostics complete; wrote {len(paths)} plots[/green]")
    return paths

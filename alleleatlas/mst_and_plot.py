"""Build MST from k-NN graph and plot components arranged on canvas.

The script loads a k-NN sparse matrix (.npz), builds a Minimum Spanning Tree
for each connected component (using the symmetric graph), computes a layout
for each component (spring layout), and draws them arranged on a larger
canvas so components don't overlap.
"""
from __future__ import annotations

import numpy as np
from typing import Optional
import math
import random
import scipy.sparse as sp
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import typer
from rich.console import Console


def read_labels_from_tsv(path: str, n_total: int, preferred: str = 'ParaC serovars'):
    """Read labels from a TSV/CSV metadata file.

    Strategy:
    - try reading with pandas (tab-separated)
    - if preferred column exists and df length == n_total -> return column values in order
    - else if df length == n_total -> pick the first non-ID column as labels
    - else if 'ID' column exists and can be parsed as integer indices -> map index->label
    - otherwise return None
    """
    try:
        df = pd.read_csv(path, sep='\t', dtype=str, low_memory=False)
    except Exception:
        try:
            df = pd.read_csv(path, dtype=str, low_memory=False)
        except Exception:
            return None

    cols = list(df.columns)
    # prefer exact preferred column
    if preferred in df.columns and len(df) == n_total:
        return df[preferred].astype(str).values

    # try common serovar-like names
    for cand in ['Serovar', 'serovar', 'serovars', 'ParaC serovars', 'ParaC_serovars', 'group', 'Group']:
        if cand in df.columns and len(df) == n_total:
            return df[cand].astype(str).values

    # if row-count matches, choose the first non-ID column
    if len(df) == n_total:
        for c in cols:
            if c.lower() not in ('id', 'index'):
                return df[c].astype(str).values
    # try mapping by integer index-like columns when df length != n_total
    index_candidates = [c for c in df.columns if c.lower() in ('index', 'id', 'sample id', 'sampleid')]
    if len(df) != n_total and len(index_candidates) > 0:
        # create labels array of length n_total filled with None
        res = [None] * n_total
        # choose label column
        lab_col = None
        for cand in (preferred, 'Serovar', 'serovar', 'serovars', 'ParaC serovars'):
            if cand in df.columns:
                lab_col = cand
                break
        if lab_col is None:
            # fallback to last column
            lab_col = df.columns[-1]
        for _, row in df.iterrows():
            for ic in index_candidates:
                try:
                    idx_val = row[ic]
                    if pd.isna(idx_val):
                        continue
                    idx = int(str(idx_val).strip())
                except Exception:
                    continue
                # support 1-based indices in TSV
                if 1 <= idx <= n_total:
                    res[idx - 1] = str(row[lab_col])
                elif 0 <= idx < n_total:
                    res[idx] = str(row[lab_col])
        return np.array(res, dtype=object)
    return None

console = Console()
app = typer.Typer()


def collapse_nodes_by_threshold(sub: sp.spmatrix, threshold: Optional[float]):
    """Collapse nodes in subgraph connected by edges with weight <= threshold.

    Returns mapping node_index -> cluster_id and number of clusters.
    """
    if threshold is None:
        return {i: i for i in range(sub.shape[0])}, sub.shape[0]
    # build graph of edges with weight <= threshold
    coo = sub.tocoo()
    rows = coo.row
    cols = coo.col
    data = coo.data
    G = nx.Graph()
    G.add_nodes_from(range(sub.shape[0]))
    for i, j, w in zip(rows, cols, data):
        if w <= threshold:
            G.add_edge(i, j)
    clusters = list(nx.connected_components(G))
    mapping = {}
    for cid, comp in enumerate(clusters):
        for v in comp:
            mapping[v] = cid
    return mapping, len(clusters)


def contract_subgraph(sub: sp.spmatrix, mapping: dict, nclusters: int) -> sp.csr_matrix:
    """Build contracted adjacency between clusters taking min inter-edge weight."""
    import collections
    edges = collections.defaultdict(lambda: np.inf)
    coo = sub.tocoo()
    for i, j, w in zip(coo.row, coo.col, coo.data):
        ci = mapping[i]
        cj = mapping[j]
        if ci == cj:
            continue
        if w < edges[(ci, cj)]:
            edges[(ci, cj)] = w
            edges[(cj, ci)] = w
    rows = []
    cols = []
    data = []
    for (i, j), w in edges.items():
        rows.append(i)
        cols.append(j)
        data.append(w)
    if len(rows) == 0:
        return sp.csr_matrix((nclusters, nclusters))
    return sp.coo_matrix((data, (rows, cols)), shape=(nclusters, nclusters)).tocsr()


def load_knn(path: str) -> sp.csr_matrix:
    return sp.load_npz(path)


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


def plot_mst(
    knn: str,
    out: str = 'mst.png',
    node_size: int = 30,
    colormap: str = 'tab10',
    layout_choice: str = 'spring',
    collapse_threshold: Optional[float] = None,
    spacing: float = 1.0,
    curved_edges: bool = False,
    spring_iterations: Optional[int] = None,
    spring_k: Optional[float] = None,
    labels_path: Optional[str] = None,
    use_embedding_init: bool = False,
    embedding_method: str = 'mds',
    embedding_mode: str = 'shortest_path',
    embedding_coords: Optional[np.ndarray] = None,
    embed_as_layout: bool = False,
    edge_prune_quantile: Optional[float] = 0.9,
):
    """Programmatic entrypoint: build MSTs, optionally collapse nodes, layout, and save image.

    Parameters
    - knn: path to scipy sparse npz
    - out: output image path
    - node_size: node drawing size
    - colormap: matplotlib colormap name
    - layout_choice: 'spring'|'circular'|'spectral'|'forceatlas2'|'precomputed'
    - collapse_threshold: merge nodes with inter-edge weight <= threshold
    - spacing: space multiplier between components
    """
    M = load_knn(knn)
    console.print(f"Loaded k-NN: shape={M.shape}, nnz={M.nnz}")
    # symmetric adjacency
    A = (M + M.T) / 2
    # find connected components on boolean adjacency
    ncomp, labels = sp.csgraph.connected_components(A, directed=False)
    console.print(f"Found {ncomp} connected components")
    comps = []
    for comp_id in range(ncomp):
        nodes = np.where(labels == comp_id)[0]
        if len(nodes) == 0:
            continue
        sub = A[nodes[:, None], nodes]
        # collapse nodes into supernodes if requested
        mapping, nclusters = collapse_nodes_by_threshold(sub, collapse_threshold)
        if nclusters < sub.shape[0]:
            C = contract_subgraph(sub, mapping, nclusters)
            try:
                Gsub = nx.from_scipy_sparse_array(C)
            except Exception:
                Gsub = nx.from_scipy_sparse_matrix(C)
        else:
            try:
                Gsub = nx.from_scipy_sparse_array(sub)
            except Exception:
                Gsub = nx.from_scipy_sparse_matrix(sub)
        if Gsub.number_of_edges() > 0:
            T = nx.minimum_spanning_tree(Gsub)
        else:
            T = Gsub
        comps.append((nodes, T, mapping, nclusters))

    # layout each component separately and compute bounding boxes
    layouts = []
    sizes = []
    # use deterministic random for small jitters/rotations
    rand = random.Random(42)
    for comp_idx, (nodes, T, mapping, nclusters) in enumerate(comps):
        # tuned spring layout: use spectral init and scale iterations by component size
        n_nodes = max(1, T.number_of_nodes())
        if layout_choice == 'spring':
            # base k heuristic unless overridden
            k = spring_k if spring_k is not None else max(0.1, math.sqrt(n_nodes) / 5.0)
            # dynamic iterations: scale with n_nodes but cap to avoid runaway
            iters = spring_iterations if spring_iterations is not None else max(200, min(2000, int(n_nodes * 2)))
            # determine position initialization: prefer embedding-based init when requested
            pos0 = None
            if use_embedding_init:
                try:
                    # lazy import to avoid dependency cycles; compute_embedding expects full indices
                    from approx.analysis.mds_from_knn import compute_embedding
                    n_total = M.shape[0]
                    # use provided embedding coords if supplied
                    if embedding_coords is not None:
                        Y_all = embedding_coords
                    else:
                        # compute global embedding for all nodes (fast for small n)
                        Y_all = compute_embedding(M, np.arange(n_total), method=embedding_method, mode=embedding_mode, random_state=42)
                    # build pos0 mapping for this component
                    if nclusters < (T.number_of_nodes()):
                        # contracted clusters -> average member coordinates
                        cluster_to_locals = {}
                        for local_idx, cid in mapping.items():
                            cluster_to_locals.setdefault(cid, []).append(local_idx)
                        pos0 = {}
                        for node in T.nodes():
                            locals_in_cluster = cluster_to_locals.get(int(node), [])
                            globals_in_cluster = [int(nodes[li]) for li in locals_in_cluster]
                            if len(globals_in_cluster) == 0:
                                pos0[node] = (0.0, 0.0)
                            else:
                                coords = Y_all[globals_in_cluster]
                                pos0[node] = (float(coords[:, 0].mean()), float(coords[:, 1].mean()))
                    else:
                        pos0 = {n: (float(Y_all[int(nodes[int(n)]), 0]), float(Y_all[int(nodes[int(n)]), 1])) for n in T.nodes()}
                except Exception:
                    pos0 = None
            # fallback to spectral init if no embedding
            if pos0 is None:
                try:
                    pos0 = nx.spectral_layout(T)
                except Exception:
                    pos0 = None
            # If requested, use embedding positions directly as the layout (no spring relax)
            if embed_as_layout and pos0 is not None:
                pos = pos0
            else:
                pos = nx.spring_layout(T, k=k, iterations=iters, pos=pos0, seed=42 + comp_idx)
        elif layout_choice == 'circular':
            pos = nx.circular_layout(T)
        elif layout_choice == 'spectral':
            pos = nx.spectral_layout(T)
        elif layout_choice == 'forceatlas2':
            # optional dependency: fa2 or other forceatlas implementation
            try:
                from fa2 import ForceAtlas2
            except Exception:
                console.print('[yellow]forceatlas2 not available, falling back to spring layout[/yellow]')
                pos = nx.spring_layout(T)
            else:
                fa = ForceAtlas2()
                pos = fa.forceatlas2_networkx_layout(T)
        else:
            pos = nx.spring_layout(T)
        xs = np.array([p[0] for p in pos.values()])
        ys = np.array([p[1] for p in pos.values()])
        # normalize to non-zero extents
        w = xs.max() - xs.min() if len(xs) > 0 else 1.0
        h = ys.max() - ys.min() if len(ys) > 0 else 1.0
        # small deterministic jitter to help break exact overlaps
        jitter_scale = max(w, h) * 1e-3
        for n in pos:
            jx = (rand.random() - 0.5) * jitter_scale
            jy = (rand.random() - 0.5) * jitter_scale
            pos[n] = (pos[n][0] + jx, pos[n][1] + jy)
        layouts.append((T, pos))
        sizes.append((w, h))

    # grid packing: arrange components in a near-square grid of cells
    ncomps = len(layouts)
    if ncomps == 0:
        console.print('[yellow]No components to plot[/yellow]')
        return
    cols = int(math.ceil(math.sqrt(ncomps)))
    rows = int(math.ceil(ncomps / cols))
    padding = 40
    # compute cell sizes proportional to component extents
    cell_w = max(300, int(1024 / cols))
    cell_h = max(300, int(1024 / rows))
    positions = []
    idx = 0
    for r in range(rows):
        for c0 in range(cols):
            if idx >= ncomps:
                break
            px = c0 * (cell_w + padding)
            py = r * (cell_h + padding)
            positions.append((px, py, cell_w, cell_h))
            idx += 1
    total_w = cols * (cell_w + padding)
    total_h = rows * (cell_h + padding)
    fig, ax = plt.subplots(figsize=(max(8, total_w / 150), max(4, total_h / 150)))
    ax.set_axis_off()
    cmap = plt.get_cmap(colormap)
    for idx, ((nodes, T, mapping, nclusters), (T2, pos), (px, py, sw, sh)) in enumerate(zip(comps, layouts, positions)):
        # transform pos into box [px, px+sw] x [py, py+sh]
        xs = np.array([p[0] for p in pos.values()])
        ys = np.array([p[1] for p in pos.values()])
        if xs.max() - xs.min() == 0:
            scale_x = 1.0
        else:
            scale_x = sw / (xs.max() - xs.min())
        if ys.max() - ys.min() == 0:
            scale_y = 1.0
        else:
            scale_y = sh / (ys.max() - ys.min())
        # optionally rotate the component slightly to reduce overlapping long edges
        angle = (idx % 5 - 2) * 0.05
        cos_a = math.cos(angle)
        sin_a = math.sin(angle)
        # center
        cx = px + sw / 2
        cy = py + sh / 2
        node_pos = {}
        for i, n in enumerate(T2.nodes()):
            ox, oy = pos[n]
            # rotate around origin of component before scaling/centering
            dx = ox - xs.mean()
            dy = oy - ys.mean()
            rx = dx * cos_a - dy * sin_a
            ry = dx * sin_a + dy * cos_a
            nxp = cx + rx * scale_x
            nyp = cy + ry * scale_y
            node_pos[n] = (nxp, nyp)
        # Optionally prune the longest MST edges to avoid long bridges that tangle clusters
        if edge_prune_quantile is not None and T2.number_of_edges() > 0:
            weights = [d.get('weight', 1.0) for _, _, d in T2.edges(data=True)]
            try:
                threshold = float(np.quantile(weights, edge_prune_quantile))
            except Exception:
                threshold = None
            if threshold is not None:
                remove = [(u, v) for u, v, d in T2.edges(data=True) if d.get('weight', 1.0) > threshold]
                if len(remove) > 0:
                    T2.remove_edges_from(remove)

        # draw
        # determine node colors: if labels provided, color by majority label per node/cluster
        color_hex = mcolors.to_hex(cmap(idx % cmap.N))
        node_colors = None
        if labels_path is not None:
            # load labels once for the whole graph
            if 'labels_all' not in locals():
                labels_all = read_labels_from_tsv(labels_path, M.shape[0])
                if labels_all is None:
                    labels_all = None
                else:
                    unique_labels = list(np.unique(labels_all))
                    label2idx = {lab: i for i, lab in enumerate(unique_labels)}
                    cmap_labels = plt.get_cmap(colormap)
            if labels_all is not None:
                node_colors = []
                # mapping maps local index in sub -> cluster id; nodes is array of global indices
                if nclusters < (T.number_of_nodes()):
                    # contracted: T2 nodes are cluster ids; find member locals per cluster
                    cluster_to_locals = {}
                    for local_idx, cid in mapping.items():
                        cluster_to_locals.setdefault(cid, []).append(local_idx)
                    for n in T2.nodes():
                        locals_in_cluster = cluster_to_locals.get(int(n), [])
                        globals_in_cluster = [int(nodes[li]) for li in locals_in_cluster]
                        # majority label
                        labs = [labels_all[g] for g in globals_in_cluster]
                        from collections import Counter

                        if len(labs) == 0:
                            node_colors.append(color_hex)
                        else:
                            maj = Counter(labs).most_common(1)[0][0]
                            node_colors.append(mcolors.to_hex(cmap_labels(label2idx[maj] % cmap_labels.N)))
                else:
                    # not contracted: T2 nodes are local indices mapping to nodes array
                    for n in T2.nodes():
                        gidx = int(nodes[int(n)])
                        lab = labels_all[gidx]
                        node_colors.append(mcolors.to_hex(cmap_labels(label2idx[lab] % cmap_labels.N)))
        # draw edges: optionally draw curved FancyArrowPatch per-edge (slower but clearer),
        # otherwise use networkx fast draw (LineCollection)
        if curved_edges:
            from matplotlib.patches import FancyArrowPatch
            # draw each edge as a curved patch; rad depends on node-pair separation
            for u, v, d in T2.edges(data=True):
                p1 = node_pos[u]
                p2 = node_pos[v]
                # compute a radius that increases with euclidean distance to make longer
                # edges slightly more curved (breaks parallel lines)
                dist = math.hypot(p2[0] - p1[0], p2[1] - p1[1])
                # normalize rad into a modest range
                rad = 0.02 + min(0.3, dist / max(sw, sh) * 0.08)
                arrow = FancyArrowPatch(p1, p2, connectionstyle=f"arc3,rad={rad}",
                                        linewidth=0.6, alpha=0.7, color=color_hex,
                                        arrowstyle='-', shrinkA=0.0, shrinkB=0.0)
                ax.add_patch(arrow)
        else:
            try:
                nx.draw_networkx_edges(T2, node_pos, ax=ax, alpha=0.6, width=0.7, connectionstyle='arc3,rad=0.08')
            except TypeError:
                # older networkx may ignore connectionstyle for certain backends
                nx.draw_networkx_edges(T2, node_pos, ax=ax, alpha=0.6, width=0.7)
        if node_colors is None:
            nx.draw_networkx_nodes(T2, node_pos, node_size=node_size, node_color=color_hex, ax=ax, linewidths=0.5)
        else:
            nx.draw_networkx_nodes(T2, node_pos, node_size=node_size, node_color=node_colors, ax=ax, linewidths=0.5)
    # if labels were provided, draw a legend mapping colors to labels
    if labels_path is not None and 'labels_all' in locals() and labels_all is not None:
        # create a small legend figure appended below main figure
        # map unique labels to colors
        # filter out None/missing labels
        labels_filtered = [l for l in labels_all if (l is not None and str(l).strip() != '')]
        if len(labels_filtered) == 0:
            unique_labels = []
        else:
            unique_labels, counts = np.unique(labels_filtered, return_counts=True)
            # sort labels by frequency descending
            order = np.argsort(counts)[::-1]
            unique_labels = list(unique_labels[order])
        label2idx = {lab: i for i, lab in enumerate(unique_labels)}
        cmap_labels = plt.get_cmap(colormap)
        # build legend patches
        from matplotlib.patches import Patch

        # compact legend if too many labels
        max_legend = 20
        legend_handles = []
        if len(unique_labels) == 0:
            legend_handles = []
        elif len(unique_labels) <= max_legend:
            for lab, i in label2idx.items():
                legend_handles.append(Patch(facecolor=mcolors.to_hex(cmap_labels(i % cmap_labels.N)), label=lab))
        else:
            # top-k + Other
            for lab in unique_labels[:max_legend - 1]:
                i = label2idx[lab]
                legend_handles.append(Patch(facecolor=mcolors.to_hex(cmap_labels(i % cmap_labels.N)), label=lab))
            legend_handles.append(Patch(facecolor='#999999', label='Other'))
        # place legend on the right
        ax_leg = fig.add_axes([0.85, 0.1, 0.12, 0.8])
        ax_leg.axis('off')
        ax_leg.legend(handles=legend_handles, loc='center')
    plt.savefig(out, dpi=150)
    console.print(f"Wrote MST plot: {out}")


@app.command()
def cli(
    knn: str = typer.Argument(..., help='k-NN .npz path'),
    out: str = typer.Argument('mst.png'),
    node_size: int = typer.Option(30, help='Node size in plot'),
    colormap: str = typer.Option('tab10', help='Matplotlib colormap name'),
    layout_choice: str = typer.Option('spring', help='Layout: spring|circular|spectral|forceatlas2'),
    collapse_threshold: Optional[float] = typer.Option(None, help='Collapse nodes with inter-edge weight <= threshold'),
    spacing: float = typer.Option(1.0, help='Component spacing multiplier'),
    spring_iterations: Optional[int] = typer.Option(None, help='Override spring layout iteration count'),
    spring_k: Optional[float] = typer.Option(None, help='Override spring layout k parameter'),
    curved_edges: bool = typer.Option(False, help='Draw curved edges (slower)') ,
):
    # When invoked programmatically (tests call cli() directly), Typer Option defaults
    # may be passed through as OptionInfo objects. Coerce to plain values in that case.
    if isinstance(node_size, typer.models.OptionInfo):
        node_size = node_size.default
    if isinstance(colormap, typer.models.OptionInfo):
        colormap = colormap.default
    if isinstance(layout_choice, typer.models.OptionInfo):
        layout_choice = layout_choice.default
    if isinstance(collapse_threshold, typer.models.OptionInfo):
        collapse_threshold = None
    if isinstance(spacing, typer.models.OptionInfo):
        spacing = spacing.default
    if isinstance(spring_iterations, typer.models.OptionInfo):
        spring_iterations = spring_iterations.default
    if isinstance(spring_k, typer.models.OptionInfo):
        spring_k = spring_k.default
    if isinstance(curved_edges, typer.models.OptionInfo):
        curved_edges = curved_edges.default
    plot_mst(
        knn,
        out,
        node_size=node_size,
        colormap=colormap,
        layout_choice=layout_choice,
        collapse_threshold=collapse_threshold,
        spacing=spacing,
        curved_edges=curved_edges,
        spring_iterations=spring_iterations,
        spring_k=spring_k,
    )


if __name__ == '__main__':
    app()

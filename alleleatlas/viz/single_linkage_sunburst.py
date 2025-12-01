"""Single-linkage sunburst plot for k-NN graphs."""

from __future__ import annotations

import math
from typing import Dict, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from matplotlib.patches import Wedge
from networkx.utils import UnionFind
from rich.console import Console

console = Console()


def load_knn(path: str) -> sp.csr_matrix:
    """Load a sparse k-NN matrix from .npz."""
    return sp.load_npz(path)


def single_linkage_labels_from_mst(A: sp.spmatrix, thresholds: Iterable[float]) -> Dict[float, np.ndarray]:
    """Return cluster labels for each threshold using a single MST sweep."""
    thresholds = sorted({float(t) for t in thresholds})
    if not thresholds:
        raise ValueError("Provide at least one threshold for sunburst clustering.")

    n = A.shape[0]
    if n == 0:
        return {t: np.array([], dtype=np.int32) for t in thresholds}

    T = sp.csgraph.minimum_spanning_tree(A.tocsr())
    coo = T.tocoo()
    edges = np.vstack([coo.data, coo.row, coo.col]).T
    if edges.size:
        edges = edges[edges[:, 0].argsort()]
    else:
        edges = edges.reshape(0, 3)

    uf = UnionFind(range(n))
    labels_per_t: Dict[float, np.ndarray] = {}
    edge_idx = 0

    for thr in thresholds:
        while edge_idx < len(edges) and edges[edge_idx, 0] <= thr:
            _, u, v = edges[edge_idx]
            uf.union(int(u), int(v))
            edge_idx += 1

        root_to_id = {}
        labels = np.empty(n, dtype=np.int32)
        next_id = 0
        for i in range(n):
            root = uf[i]
            cid = root_to_id.get(root)
            if cid is None:
                cid = next_id
                root_to_id[root] = cid
                next_id += 1
            labels[i] = cid
        labels_per_t[thr] = labels

    return labels_per_t


def plot_single_linkage_sunburst(
    knn: str,
    thresholds: Iterable[float],
    out: str = "cluster_sunburst.png",
    center_radius: float = 0.2,
    cmap_name: str = "tab20",
    min_fraction: float = 0.0,
) -> None:
    """Draw a sunburst where each ring is a single-linkage cut."""
    thresholds = list(thresholds)
    if not thresholds:
        raise ValueError("Provide at least one threshold for the sunburst plot.")
    if not (0.0 <= min_fraction < 1.0):
        raise ValueError("min_fraction must be in [0, 1).")

    console.print(f"[blue]Loading k-NN graph for sunburst from:[/blue] {knn}")
    M = load_knn(knn)
    A = (M + M.T) / 2

    labels_per_t = single_linkage_labels_from_mst(A, thresholds)
    n = A.shape[0]
    console.print(f"[blue]Profiles (nodes) in sunburst:[/blue] {n}")

    # Rings: inner = coarsest (largest threshold), outer = finest (smallest threshold)
    levels = sorted(labels_per_t.keys(), reverse=True)
    labels_levels = [labels_per_t[t] for t in levels]
    n_levels = len(levels)

    root = {"name": "all", "count": n, "children": {}, "is_other": False}
    for idx in range(n):
        node = root
        for lvl, lbls in enumerate(labels_levels):
            cid = int(lbls[idx])
            key = (lvl, cid)
            child = node["children"].get(key)
            if child is None:
                child = {
                    "name": f"t={levels[lvl]:g}, c={cid}",
                    "count": 0,
                    "children": {},
                    "is_other": False,
                }
                node["children"][key] = child
            child["count"] += 1
            node = child

    def collapse_small(node, level: int) -> None:
        children_items = list(node["children"].items())
        if not children_items:
            return
        total = sum(ch["count"] for _, ch in children_items)
        if total <= 0:
            return

        keep = []
        other_count = 0
        for key, ch in children_items:
            frac = ch["count"] / total if total > 0 else 0.0
            if min_fraction > 0.0 and frac < min_fraction:
                other_count += ch["count"]
            else:
                keep.append((key, ch))
        new_children = {k: v for k, v in keep}
        if other_count > 0:
            other_key = (level, f"other-{level}")
            new_children[other_key] = {
                "name": "other",
                "count": other_count,
                "children": {},
                "is_other": True,
            }
        node["children"] = new_children
        for child in node["children"].values():
            collapse_small(child, level + 1)

    if min_fraction > 0.0:
        collapse_small(root, 0)

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={"aspect": "equal"})
    cmap = plt.get_cmap(cmap_name)
    ring_width = (1.0 - center_radius) / n_levels if n_levels > 0 else 0.0

    def lighten_color(color, amount=0.5):
        r, g, b, a = mcolors.to_rgba(color)
        r = 1 - (1 - r) * (1 - amount)
        g = 1 - (1 - g) * (1 - amount)
        b = 1 - (1 - b) * (1 - amount)
        return (r, g, b, a)

    if n > 0:
        center_patch = Wedge(
            (0.0, 0.0),
            r=center_radius,
            theta1=0.0,
            theta2=360.0,
            width=center_radius,
            facecolor="#dddddd",
            edgecolor="white",
            linewidth=0.5,
        )
        ax.add_patch(center_patch)
        ax.text(0, 0, f"N={n}", ha="center", va="center", fontsize=9)

    def draw_children(node, level, theta_start, theta_end, parent_color=None):
        children = sorted(
            node["children"].values(),
            key=lambda ch: (ch.get("is_other", False), -ch["count"]),
        )
        if not children or level >= n_levels:
            return
        total = sum(ch["count"] for ch in children)
        if total <= 0:
            return

        inner_r = center_radius + level * ring_width
        outer_r = inner_r + ring_width
        angle = theta_start

        for idx, ch in enumerate(children):
            frac = ch["count"] / total
            th0 = angle
            th1 = angle + frac * (theta_end - theta_start)

            if ch.get("is_other"):
                color = "#d3d3d3"
                base_color = parent_color
            else:
                if level == 0:
                    base_color = cmap(idx / max(1, len(children)))
                else:
                    base_color = parent_color if parent_color is not None else cmap(0.0)
                mix = 0.3 + 0.5 * (level / max(1, n_levels - 1))
                color = lighten_color(base_color, amount=mix)

            wedge = Wedge(
                (0.0, 0.0),
                r=outer_r,
                theta1=math.degrees(th0),
                theta2=math.degrees(th1),
                width=ring_width * 0.98 if ring_width > 0 else None,
                facecolor=color,
                edgecolor="white",
                linewidth=0.5,
            )
            ax.add_patch(wedge)
            draw_children(
                ch,
                level + 1,
                th0,
                th1,
                parent_color=base_color,
            )
            angle = th1

    draw_children(root, level=0, theta_start=0.0, theta_end=2 * math.pi, parent_color=None)

    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    console.print(f"[green]Wrote single-linkage sunburst plot:[/green] {out}")

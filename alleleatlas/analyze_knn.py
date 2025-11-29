"""Analyze k-NN graphs: distance distributions, connectivity, and plots.

Produces simple matplotlib plots showing distance histogram, degree
distribution, and number of connected components as a function of a
distance cutoff applied to the k-NN edges.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import typer
from rich.console import Console
import matplotlib.pyplot as plt
import json
import math

console = Console()
app = typer.Typer()


def load_knn(npz_path: str) -> sp.csr_matrix:
    return sp.load_npz(npz_path)


def edge_distances(M: sp.spmatrix) -> np.ndarray:
    # return distances for non-zero entries (upper triangle)
    coo = M.tocoo()
    mask = coo.row < coo.col
    return coo.data[mask]


def _find_two_elbows(distances: np.ndarray) -> tuple[float, float]:
    """Simple two-elbow finder using a Kneedle-like heuristic on sorted distances."""
    if len(distances) == 0:
        return (0.0, 0.0)
    d = np.sort(distances)

    def _knee(vals: np.ndarray, offset: int = 0) -> int:
        if len(vals) == 0:
            return offset
        x = np.linspace(0, 1, len(vals))
        spread = np.ptp(vals) if len(vals) > 0 else 0.0
        y = (vals - vals.min()) / (spread if spread > 0 else 1.0)
        scores = y - x  # convex elbow
        return offset + int(np.argmax(scores))

    idx1 = _knee(d, 0)
    idx2 = _knee(d[idx1 + 1 :], idx1 + 1) if idx1 + 1 < len(d) else idx1
    fine = float(d[idx1])
    coarse = float(d[idx2])
    # fallback: if elbows collapse or sit at the max, use quantiles instead
    if fine == coarse or coarse >= d[-1] or fine >= d[-1]:
        q50 = float(np.quantile(d, 0.5))
        q80 = float(np.quantile(d, 0.8))
        fine, coarse = q50, q80
    return fine, coarse


def _component_curve(A: sp.spmatrix, cutoffs: np.ndarray):
    """Compute #components and largest component size for each cutoff."""
    coo = A.tocoo()
    comps = []
    giant = []
    for c in cutoffs:
        mask = coo.data <= c
        if not np.any(mask):
            comps.append(A.shape[0])
            giant.append(1)
            continue
        G = sp.coo_matrix((coo.data[mask], (coo.row[mask], coo.col[mask])), shape=A.shape)
        ncomp, labels = sp.csgraph.connected_components(G, directed=False)
        comps.append(int(ncomp))
        if labels is not None and len(labels) > 0:
            counts = np.bincount(labels)
            giant.append(int(counts.max()))
        else:
            giant.append(1)
    return np.array(comps), np.array(giant)


def _elbows_from_curve(values: np.ndarray, cutoffs: np.ndarray) -> tuple[float, float]:
    """Find two knees on a decreasing curve (components vs cutoff)."""
    if len(values) == 0 or len(values) != len(cutoffs):
        return (cutoffs[0] if len(cutoffs) else 0.0, cutoffs[-1] if len(cutoffs) else 0.0)
    # normalize
    v = values.astype(float)
    if v.max() == v.min():
        return (cutoffs[0], cutoffs[-1])
    y = (v - v.min()) / (v.max() - v.min())
    x = np.linspace(0, 1, len(v))
    # For decreasing curves, use (1 - y) to make an increasing convex-like shape
    y2 = 1 - y
    scores = y2 - x
    idx1 = int(np.argmax(scores))
    idx2 = idx1 + int(np.argmax(scores[idx1 + 1:])) + 1 if idx1 + 1 < len(scores) else idx1
    return float(cutoffs[idx1]), float(cutoffs[idx2])


@app.command()
def cli(
    knn: str = typer.Argument(..., help='k-NN .npz file (scipy sparse)') ,
    out_prefix: str = typer.Argument('knn_analysis', help='Output prefix for plots') ,
    cutoff_steps: int = 50,
):
    run_analysis(knn, out_prefix, cutoff_steps)


def run_analysis(knn: str, out_prefix: str = 'knn_analysis', cutoff_steps: int = 50):
    """Programmatic entrypoint for analysis (callable from tests).

    The Typer CLI ``cli`` wraps this function.
    """
    M = load_knn(knn)
    console.print(f"Loaded k-NN: shape={M.shape}, nnz={M.nnz}")
    dists = edge_distances(M)
    # distance-based elbows (fallback)
    elbow_dist_fine, elbow_dist_coarse = _find_two_elbows(dists)
    # component curve
    A = (M + M.T) / 2
    min_d, max_d = (float(dists.min()) if len(dists) else 0.0, float(dists.max()) if len(dists) else 1.0)
    # choose cutoffs spanning quantiles to avoid extreme tails
    lo = float(np.quantile(dists, 0.05)) if len(dists) else min_d
    hi = float(np.quantile(dists, 0.95)) if len(dists) else max_d
    cutoffs = np.linspace(lo, hi, cutoff_steps)
    comp_counts, giant_sizes = _component_curve(A, cutoffs)
    elbow_comp_fine, elbow_comp_coarse = _elbows_from_curve(comp_counts, cutoffs)
    # choose component-based elbows if they differ; otherwise fallback to distance-based
    elbow1 = elbow_comp_fine if not math.isclose(elbow_comp_fine, elbow_comp_coarse) else elbow_dist_fine
    elbow2 = elbow_comp_coarse if not math.isclose(elbow_comp_fine, elbow_comp_coarse) else elbow_dist_coarse
    console.print(f"Elbow thresholds (approx): coarse={elbow2:.4f}, fine={elbow1:.4f}")
    with open(f"{out_prefix}_elbows.json", "w") as fh:
        json.dump(
            {
                "fine": elbow1,
                "coarse": elbow2,
                "distance_elbows": {"fine": elbow_dist_fine, "coarse": elbow_dist_coarse},
                "component_elbows": {"fine": elbow_comp_fine, "coarse": elbow_comp_coarse},
            },
            fh,
            indent=2,
        )
    console.print(f"Wrote {out_prefix}_elbows.json")
    # plot histogram of distances
    plt.figure()
    plt.hist(dists, bins=100)
    plt.title('k-NN edge distance distribution')
    plt.xlabel('distance')
    plt.ylabel('count')
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_dist_hist.png")
    console.print(f"Wrote {out_prefix}_dist_hist.png")

    # degree distribution
    deg = np.array(M.sum(axis=1)).reshape(-1)
    plt.figure()
    plt.hist(deg, bins=50)
    plt.title('Node degree distribution (sum of weights)')
    plt.xlabel('degree')
    plt.ylabel('count')
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_deg_hist.png")
    console.print(f"Wrote {out_prefix}_deg_hist.png")

    # components vs cutoff plot
    plt.figure()
    plt.plot(cutoffs, comp_counts, marker='o')
    plt.axvline(elbow1, color='orange', linestyle='--', label='fine')
    plt.axvline(elbow2, color='red', linestyle='--', label='coarse')
    plt.title('Connected components vs cutoff')
    plt.xlabel('distance cutoff')
    plt.ylabel('# components')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_components.png")
    console.print(f"Wrote {out_prefix}_components.png")
    return elbow1, elbow2

if __name__ == '__main__':
    app()

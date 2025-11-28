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

console = Console()
app = typer.Typer()


def load_knn(npz_path: str) -> sp.csr_matrix:
    return sp.load_npz(npz_path)


def edge_distances(M: sp.spmatrix) -> np.ndarray:
    # return distances for non-zero entries (upper triangle)
    coo = M.tocoo()
    mask = coo.row < coo.col
    return coo.data[mask]


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

    # number of connected components for a range of cutoffs
    from scipy.sparse.csgraph import connected_components
    if dists.size == 0:
        console.print('[yellow]No edges found in k-NN; skipping component sweep[/yellow]')
        return
    cutoffs = np.linspace(dists.min(), dists.max(), cutoff_steps)
    comps = []
    for c in cutoffs:
        # keep edges with distance <= c
        Mcut = M.copy()
        Mcut.data = (Mcut.data <= c).astype(int)
        ncomp, labels = connected_components(Mcut, directed=False)
        comps.append(ncomp)
    plt.figure()
    plt.plot(cutoffs, comps)
    plt.xlabel('distance cutoff')
    plt.ylabel('n_components')
    plt.title('Components vs distance cutoff')
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_components.png")
    console.print(f"Wrote {out_prefix}_components.png")


if __name__ == '__main__':
    app()

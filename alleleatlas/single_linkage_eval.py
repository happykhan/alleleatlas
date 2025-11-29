"""Single-linkage clustering evaluation on k-NN graphs.

We approximate legacy HierCC-style plots by:
- building an MST on the symmetric k-NN graph
- sweeping distance thresholds (quantiles, fixed step, or explicit list)
- recording cluster labels, component counts, NMI between levels, and optional silhouette
- writing a TSV summary plus a heatmap/line plot

Silhouette uses a subsampled shortest-path distance matrix to stay scalable.
"""
from __future__ import annotations

import math
from typing import List, Tuple, Optional
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy.sparse as sp
import typer
from sklearn.metrics import normalized_mutual_info_score
from sklearn.metrics import silhouette_score
from rich.console import Console
from networkx.utils import UnionFind

console = Console()
app = typer.Typer()


def load_knn(path: str) -> sp.csr_matrix:
    return sp.load_npz(path)


def mst_edges(A: sp.spmatrix) -> Tuple[int, np.ndarray]:
    """Return node count and sorted (weight, u, v) edges of the MST."""
    from scipy.sparse.csgraph import minimum_spanning_tree

    T = minimum_spanning_tree(A)
    coo = T.tocoo()
    edges = np.vstack([coo.data, coo.row, coo.col]).T
    edges = edges[edges[:, 0].argsort()]
    n = A.shape[0]
    return n, edges


def thresholds_from_edges(edges: np.ndarray, n_levels: int, qmin: float, qmax: float) -> np.ndarray:
    weights = edges[:, 0]
    qs = np.linspace(qmin, qmax, n_levels)
    thr = np.quantile(weights, qs)
    thr = np.unique(thr)
    return thr


def build_thresholds(
    edges: np.ndarray,
    n_levels: int,
    qmin: float,
    qmax: float,
    fixed_thresholds: Optional[List[float]] = None,
    step: Optional[float] = None,
    tmin: Optional[float] = None,
    tmax: Optional[float] = None,
) -> np.ndarray:
    """Choose thresholds either from explicit list, fixed step, or quantiles."""
    weights = edges[:, 0]
    if fixed_thresholds:
        thr = np.unique(np.sort(np.array(fixed_thresholds, dtype=float)))
    elif step is not None:
        lo = float(weights.min()) if tmin is None else float(tmin)
        hi = float(weights.max()) if tmax is None else float(tmax)
        thr = np.arange(lo, hi + 1e-9, step, dtype=float)
    else:
        thr = thresholds_from_edges(edges, n_levels=n_levels, qmin=qmin, qmax=qmax)
    return thr


def labels_from_mst(n: int, edges: np.ndarray, thresholds: np.ndarray) -> List[np.ndarray]:
    """Single-linkage cluster labels for each threshold using a growing union-find."""
    uf = UnionFind(range(n))
    labels_list = []
    edge_idx = 0
    for thr in thresholds:
        while edge_idx < len(edges) and edges[edge_idx, 0] <= thr:
            _, u, v = edges[edge_idx]
            uf.union(int(u), int(v))
            edge_idx += 1
        # map roots to compact ints
        root_to_id = {}
        lbl = np.empty(n, dtype=np.int32)
        next_id = 0
        for i in range(n):
            r = uf[i]
            if r not in root_to_id:
                root_to_id[r] = next_id
                next_id += 1
            lbl[i] = root_to_id[r]
        labels_list.append(lbl)
    return labels_list


def nmi_matrix(labels_list: List[np.ndarray]) -> np.ndarray:
    k = len(labels_list)
    sim = np.ones((k, k), dtype=float)
    for i in range(k):
        for j in range(i + 1, k):
            sim[i, j] = normalized_mutual_info_score(labels_list[i], labels_list[j])
            sim[j, i] = sim[i, j]
    return sim


def component_stats(labels_list: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    comps = []
    giants = []
    for lab in labels_list:
        uniq, counts = np.unique(lab, return_counts=True)
        comps.append(len(uniq))
        giants.append(int(counts.max()))
    return np.array(comps), np.array(giants)


def silhouette_from_graph(
    A: sp.spmatrix,
    labels_list: List[np.ndarray],
    sample_size: int = 2000,
    random_state: int = 42,
) -> np.ndarray:
    """Approximate silhouette using shortest-path distances on a subsample.

    Computing full pairwise distances for 70k+ nodes is infeasible; we instead:
    - pick a subsample of nodes
    - run Dijkstra from each sampled node to all nodes, then subset to samples
    - compute silhouette_score on the sampled distance matrix for each level
    """
    n = A.shape[0]
    sample_size = min(sample_size, n)
    rng = np.random.default_rng(random_state)
    sample_idx = np.sort(rng.choice(n, size=sample_size, replace=False))
    console.print(f"[blue]Computing silhouette on subsample of {sample_size} nodes[/blue]")
    dist = sp.csgraph.dijkstra(A, directed=False, indices=sample_idx)
    dist = dist[:, sample_idx]
    sil = []
    for lab in labels_list:
        labs = lab[sample_idx]
        uniq = np.unique(labs)
        if 2 <= uniq.size < labs.size:
            try:
                s = silhouette_score(dist, labs, metric="precomputed")
            except Exception:
                s = 0.0
        else:
            s = 0.0
        sil.append(float(s))
    return np.array(sil, dtype=float)


def infer_loci_count(profile_path: Optional[str]) -> Optional[int]:
    """Best-effort allele/locus count from a profile/parquet path."""
    if profile_path is None:
        return None
    try:
        p = Path(profile_path)
        if p.suffix == ".parquet":
            import pandas as pd

            df = pd.read_parquet(p)
            if "ST" in df.columns:
                return max(0, df.shape[1] - 1)
            return max(0, df.shape[1] - 1)
        else:
            import gzip

            opener = gzip.open if str(p).endswith(".gz") else open
            with opener(p, "rt") as fh:
                for line in fh:
                    if not line:
                        continue
                    if line.startswith("#"):
                        header = line.strip().split("\t")
                        allele_cols = [i for i, h in enumerate(header) if i > 0 and not h.startswith("#")]
                        return len(allele_cols)
                    parts = line.strip().split("\t")
                    if len(parts) > 1:
                        return max(0, len(parts) - 1)
    except Exception as e:  # pragma: no cover - best effort logging
        console.print(f"[yellow]Could not infer number of loci from {profile_path}: {e}[/yellow]")
    return None


def plot_summary(
    out_prefix: str,
    thresholds: np.ndarray,
    comps: np.ndarray,
    nmi: np.ndarray,
    silhouette: Optional[np.ndarray] = None,
):
    fig, axs = plt.subplots(
        2, 2, figsize=(8, 12), gridspec_kw={"width_ratios": (12, 1), "height_ratios": (65, 35)}
    )
    heat = axs[0, 0].imshow(
        10 * (np.log10(1 - nmi + 1e-9)),
        norm=colors.TwoSlopeNorm(vmin=-30.0, vcenter=-10.0, vmax=0),
        cmap="RdBu",
        extent=[thresholds[0], thresholds[-1], thresholds[-1], thresholds[0]],
    )
    cb = fig.colorbar(heat, cax=axs[0, 1])
    line1, = axs[1, 0].plot(thresholds, comps, marker="o", label="# components")
    axs[1, 0].set_xlim([thresholds[0], thresholds[-1]])
    if silhouette is not None:
        ax2 = axs[1, 0].twinx()
        ax2.plot(thresholds, silhouette, color="tab:orange", label="silhouette")
        ax2.set_ylabel("Silhouette (subsampled)")
        # build a combined legend
        lines = [line1, ax2.lines[0]]
        labels = [ln.get_label() for ln in lines]
        axs[1, 0].legend(lines, labels, loc="upper right")
    axs[1, 1].remove()
    axs[0, 0].set_ylabel("Thresholds")
    axs[0, 0].set_xlabel("Thresholds")
    axs[1, 0].set_ylabel("# components")
    axs[1, 0].set_xlabel("Thresholds")
    cb.set_label("Normalized Mutual Information")
    cb.set_ticks([-30, -23.01, -20, -13.01, -10, -3.01, 0])
    cb.ax.set_yticklabels([">=.999", ".995", ".99", ".95", ".9", ".5", ".0"])
    plt.tight_layout()
    plt.savefig(f"{out_prefix}.pdf")
    plt.close(fig)


def run_single_linkage_eval(
    knn: str,
    out_prefix: str = "single_linkage_eval",
    n_levels: int = 40,
    qmin: float = 0.05,
    qmax: float = 0.95,
    silhouette: bool = False,
    silhouette_sample: int = 2000,
    silhouette_seed: int = 42,
    threshold: Optional[List[float]] = None,
    threshold_step: Optional[float] = None,
    threshold_min: Optional[float] = None,
    threshold_max: Optional[float] = None,
    profiles_path: Optional[str] = None,
):
    A = load_knn(knn)
    A = (A + A.T) / 2
    console.print(f"Loaded k-NN for single-linkage eval: shape={A.shape}, nnz={A.nnz}")
    loci_count = infer_loci_count(profiles_path)
    if loci_count is not None:
        console.print(f"[blue]Inferred loci/alleles:[/blue] {loci_count}")
    n, edges = mst_edges(A)
    thresholds = build_thresholds(
        edges,
        n_levels=n_levels,
        qmin=qmin,
        qmax=qmax,
        fixed_thresholds=threshold,
        step=threshold_step,
        tmin=threshold_min,
        tmax=threshold_max,
    )
    console.print(f"Evaluating {len(thresholds)} thresholds: min={thresholds.min():.4f}, max={thresholds.max():.4f}")
    labels_list = labels_from_mst(n, edges, thresholds)
    comps, giants = component_stats(labels_list)
    sim = nmi_matrix(labels_list)
    sil = None
    if silhouette:
        try:
            sil = silhouette_from_graph(A, labels_list, sample_size=silhouette_sample, random_state=silhouette_seed)
        except Exception as e:
            console.print(f"[red]Silhouette computation failed:[/red] {e}")
            sil = None

    # write TSV summary
    with open(f"{out_prefix}.tsv", "w") as fh:
        header = ["threshold", "n_components", "giant_component"]
        if sil is not None:
            header.append("silhouette_subsample")
        fh.write("\t".join(header) + "\n")
        for idx, (t, c, g) in enumerate(zip(thresholds, comps, giants)):
            row = [f"{t}", f"{c}", f"{g}"]
            if sil is not None:
                row.append(f"{sil[idx]:.5f}")
            fh.write("\t".join(row) + "\n")
    console.print(f"[green]Wrote single-linkage summary:[/green] {out_prefix}.tsv")

    plot_summary(out_prefix, thresholds, comps, sim, silhouette=sil)
    console.print(f"[green]Wrote single-linkage heatmap/line plot:[/green] {out_prefix}.pdf")


@app.command()
def cli(
    knn: str = typer.Argument(..., help="k-NN .npz path (scipy sparse)"),
    out_prefix: str = typer.Argument("single_linkage_eval"),
    n_levels: int = typer.Option(40, help="Number of thresholds to evaluate"),
    qmin: float = typer.Option(0.05, help="Lower quantile of MST edges for threshold grid"),
    qmax: float = typer.Option(0.95, help="Upper quantile of MST edges for threshold grid"),
    silhouette: bool = typer.Option(False, help="Compute silhouette on a subsample (costly but informative)"),
    silhouette_sample: int = typer.Option(2000, help="Subsample size for silhouette"),
    silhouette_seed: int = typer.Option(42, help="Random seed for silhouette subsample"),
    threshold: List[float] = typer.Option([], "--threshold", "-t", help="Explicit distance thresholds (can repeat)"),
    threshold_step: Optional[float] = typer.Option(None, help="Fixed step between thresholds (distance units)"),
    threshold_min: Optional[float] = typer.Option(None, help="Minimum threshold when using --threshold-step"),
    threshold_max: Optional[float] = typer.Option(None, help="Maximum threshold when using --threshold-step"),
    profiles_path: Optional[str] = typer.Option(None, help="Optional profile/parquet to report #alleles"),
):
    run_single_linkage_eval(
        knn,
        out_prefix,
        n_levels,
        qmin,
        qmax,
        silhouette=silhouette,
        silhouette_sample=silhouette_sample,
        silhouette_seed=silhouette_seed,
        threshold=threshold if threshold else None,
        threshold_step=threshold_step,
        threshold_min=threshold_min,
        threshold_max=threshold_max,
        profiles_path=profiles_path,
    )


if __name__ == "__main__":
    app()

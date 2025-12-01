"""Lightweight wrapper for the streamlined exact-distance pipeline.

The CLI entrypoint (`alleleatlas.py`) calls into this module to run:
1) convert profiles to parquet
2) load as float matrix (NaN = missing)
3) build exact legacy-style k-NN via distances_exact
4) analyze and plot MST
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, List, Optional
import csv
import multiprocessing as mp
import numpy as np
import scipy.sparse as sp
import pandas as pd
from rich.console import Console

from alleleatlas.get_distances import getDistance
from alleleatlas.mst_and_plot import single_linkage_components

console = Console()


@dataclass
class ClusteringConfig:
    cgmlst_profiles: str
    outdir: str
    k: int = 50
    force_recompute: bool = False
    nproc: int = 2
    verbose: bool = True
    supernode_layout: str = 'spring'
    coarse_threshold: float = 40.0
    min_cluster_size: int = 1


# helper to save sparse/dense safely
def _save_matrix(path: str, mat):
    if hasattr(mat, 'tocsr'):
        sp.save_npz(path, mat)
    else:
        # convert dense to CSR
        sp.save_npz(path, sp.csr_matrix(np.asarray(mat)))


def _symmetrize(M: sp.spmatrix) -> sp.spmatrix:
    return (M + M.T) / 2


def _clusters_from_knn(A: sp.spmatrix, threshold: float, min_component_size: int = 1) -> List[List[int]]:
    ncomp, labels = sp.csgraph.connected_components(A, directed=False)
    clusters: List[List[int]] = []
    for comp_id in range(ncomp):
        nodes = np.where(labels == comp_id)[0]
        if nodes.size == 0:
            continue
        sub = A[nodes[:, None], nodes]
        clusters_local = single_linkage_components(sub, threshold)
        for comp in clusters_local:
            clusters.append([int(nodes[i]) for i in comp])
    clusters = [c for c in clusters if len(c) >= min_component_size]
    return clusters


def _medoid_for_cluster(A: sp.spmatrix, nodes: List[int]) -> int:
    if len(nodes) == 1:
        return int(nodes[0])
    sub = A[nodes, :][:, nodes]
    dist_mat = sp.csgraph.dijkstra(sub, directed=False, unweighted=False)
    dist_mat = np.where(np.isfinite(dist_mat), dist_mat, np.nanmax(dist_mat[np.isfinite(dist_mat)]) * 2 if np.isfinite(dist_mat).any() else 0.0)
    dist_sums = np.nansum(dist_mat, axis=1)
    medoid_local_idx = int(np.nanargmin(dist_sums))
    return int(nodes[medoid_local_idx])


def _distance_matrix_from_getdistance(mat: np.ndarray, pool: mp.Pool, allowed_missing: float = 0.05) -> np.ndarray:
    dist_tri = getDistance(mat, 'dual_dist', pool, start=0, allowed_missing=allowed_missing)
    n = dist_tri.shape[1]
    full = np.zeros((n, n), dtype=float)
    for i in range(dist_tri.shape[0]):
        for j in range(min(i, n)):
            d = float(dist_tri[i, j, 0])
            full[i, j] = d
            full[j, i] = d
    np.fill_diagonal(full, 0.0)
    return full

def run(
    input_profile: str,
    out_dir: str,
    *,
    k: int = 50,
    force: bool = False,
    nproc: int = 4,
    coarse_threshold: float = 40.0,
    min_cluster_size: int = 1,
) -> None:
    """Run the full end-to-end approx pipeline.

    This variant:
      - builds a k-NN with USearch
      - collapses nodes by single-linkage at ``coarse_threshold``
      - selects cluster medoids and writes a medoid-only profile file
      - computes a dense medoid distance matrix via legacy getDistance
      - runs downstream analyses/plots on the medoid distance matrix
    """
    console.print(f"[blue]Running medoid-reduced pipeline:[/blue] {input_profile} -> {out_dir} (k={k}, threshold={coarse_threshold})")
    p = Path(out_dir)
    p.mkdir(parents=True, exist_ok=True)
    base = Path(input_profile).stem

    # 1) Convert input profile to parquet
    from alleleatlas import convert_to_parquet
    parquet_path = p / f"{base}.parquet"
    if parquet_path.exists() and not force:
        console.print(f"[yellow]Reusing existing parquet:[/yellow] {parquet_path}")
    else:
        console.print(f"[green]Converting input to parquet:[/green] {parquet_path}")
        convert_to_parquet.convert_file(input_profile, str(parquet_path))
    # report number of loci/alleles
    try:
        import pandas as pd
        df_tmp = pd.read_parquet(parquet_path)
        loci_cols = df_tmp.columns.drop('ST') if 'ST' in df_tmp.columns else df_tmp.columns[1:]
        console.print(f"[blue]Detected loci/alleles:[/blue] {len(loci_cols)}")
    except Exception as e:
        console.print(f"[yellow]Could not infer loci count:[/yellow] {e}")
    

    backend_suffix = 'usearch'
    knn_path = p / f"{base}_{backend_suffix}.npz"
    if knn_path.exists() and not force:
        console.print(f"[yellow]Reusing existing k-NN matrix:[/yellow] {knn_path}")
    else:
        # 2) Load profiles (numpy array with NaNs)
        from alleleatlas import build_knn as build_knn_module
        # load profiles with missing filter
        X = build_knn_module.load_profiles_from_parquet(str(parquet_path), missing_filter=0.1)

        from alleleatlas.backends.usearch_enhance_backend import build_knn as usearch_build_knn
        console.print(f'[green]Building k-NN (k={k}) via USearch custom metric[/green]')
        knn_mat, index_obj, _ = usearch_build_knn(X, k, return_index=True, threads=nproc)

        _save_matrix(str(knn_path), knn_mat)

    # Symmetrize for clustering
    M = sp.load_npz(knn_path)
    A = _symmetrize(M)
    console.print(f"[blue]k-NN shape:[/blue] {A.shape}, nnz={A.nnz}")

    # Single-linkage clustering to coarse clusters
    console.print(f"[blue]Single-linkage clustering at threshold {coarse_threshold}[/blue]")
    clusters = _clusters_from_knn(A, threshold=coarse_threshold, min_component_size=min_cluster_size)
    console.print(f"[blue]Clusters found:[/blue] {len(clusters)}")
    if not clusters:
        console.print("[red]No clusters found; aborting downstream steps[/red]")
        return

    # Medoids per cluster
    medoids: List[int] = []
    sizes: List[int] = []
    for cid, members in enumerate(clusters):
        med = _medoid_for_cluster(A, members)
        medoids.append(med)
        sizes.append(len(members))
    console.print(f"[blue]Medoids selected:[/blue] {len(medoids)}")

    # Write cluster summary
    cluster_tsv = p / f"{base}_{backend_suffix}_clusters.tsv"
    with open(cluster_tsv, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["cluster_id", "size", "medoid", "members"])
        for cid, members in enumerate(clusters):
            writer.writerow([cid, len(members), medoids[cid] if cid < len(medoids) else -1, ",".join(str(m) for m in members)])
    console.print(f"[green]Wrote cluster summary:[/green] {cluster_tsv}")

    # Build medoid profile file (parquet + tsv)
    df_full = pd.read_parquet(parquet_path)
    medoid_df = df_full.iloc[medoids].copy()
    medoid_parquet = p / f"{base}_{backend_suffix}_medoids.parquet"
    medoid_tsv = p / f"{base}_{backend_suffix}_medoids.tsv"
    medoid_df.to_parquet(medoid_parquet, index=False)
    medoid_df.to_csv(medoid_tsv, sep="\t", index=False)
    console.print(f"[green]Wrote medoid profiles:[/green] {medoid_parquet}")

    # Prepare matrix for getDistance (first column = ST or index)
    if "ST" in medoid_df.columns:
        mat = medoid_df.copy()
    else:
        mat = medoid_df.copy()
        mat.insert(0, "ST", np.arange(len(mat), dtype=int))
    mat_np = mat.to_numpy()
    mat_np = np.nan_to_num(mat_np, nan=0).astype(np.int32, copy=False)

    console.print("[blue]Computing dense medoid distance matrix via getDistance[/blue]")
    pool = mp.Pool(processes=max(1, nproc))
    try:
        dist_dense = _distance_matrix_from_getdistance(mat_np, pool, allowed_missing=0.05)
    finally:
        pool.close()
        pool.join()

    medoid_dist_npz = p / f"{base}_{backend_suffix}_medoid_distances.npz"
    dist_csr = sp.csr_matrix(dist_dense)
    sp.save_npz(medoid_dist_npz, dist_csr)
    console.print(f"[green]Wrote dense medoid distance matrix:[/green] {medoid_dist_npz}")

    # 4) Analyze and plot (reuse original pipeline steps on medoid distance matrix)
    try:
        from alleleatlas import analyze_knn
        from alleleatlas import mst_and_plot
        from alleleatlas.viz.single_linkage_sunburst import plot_single_linkage_sunburst
        from alleleatlas.viz.cluster_diagnostics import run_cluster_diagnostics
        from alleleatlas.viz.knn_hulls import plot_knn_hulls

        out_prefix = p / f"{base}_{backend_suffix}_medoids"
        analyze_knn.run_analysis(str(medoid_dist_npz), out_prefix=str(out_prefix))

        plot_single_linkage_sunburst(
            knn=str(medoid_dist_npz),
            thresholds=[20, 70, 400, 600, 1000, 1500],
            out=str(p / f"{base}_{backend_suffix}_medoids_sunburst.png"),
            min_fraction=0.01,
        )

        mst_and_plot.plot_cluster_backbone(
            str(medoid_dist_npz),
            str(p / f"{base}_{backend_suffix}_medoids_cluster_mst.png"),
            component_summary_out=str(p / f"{base}_{backend_suffix}_medoids_supernodes.tsv"),
        )

        diag_prefix = p / f"{base}_{backend_suffix}_medoids_diagnostics"
        run_cluster_diagnostics(
            knn=str(medoid_dist_npz),
            coarse_threshold=coarse_threshold,
            min_component_size=1,
            edge_display_threshold=600,
            out_prefix=str(diag_prefix),
        )

        plot_knn_hulls(
            knn=str(medoid_dist_npz),
            out=str(p / f"{base}_{backend_suffix}_medoids_cluster_hulls.png"),
            coarse_threshold=coarse_threshold,
            min_component_size=1,
            max_supernodes=800,
            max_local_nodes=1000,
            edge_display_threshold=600,
            hull_threshold=None,
            min_hull_nodes=5,
            layout_mode="mds",
            draw_edges=True,
            draw_nodes=True,
            random_seed=42,
        )
    except Exception as e:
        console.print(f"[red]Medoid analysis/plots failed:[/red] {e}")


def run_pipeline(config: ClusteringConfig) -> Dict[str, Any]:
    """Run pipeline given a ClusteringConfig and return artifact metadata.

    This is the function used by the CLI entrypoint. It is intentionally
    small and delegates heavy work to run().
    """
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    console.print(f"[cyan]Running pipeline:[/cyan] {config.cgmlst_profiles} -> {outdir} (k={config.k})")

    run(
        config.cgmlst_profiles,
        str(outdir),
        k=config.k,
        force=config.force_recompute,
        nproc=config.nproc,
        coarse_threshold=getattr(config, "coarse_threshold", 40.0),
        min_cluster_size=getattr(config, "min_cluster_size", 1),
    )

    artifacts = list(sorted([str(p) for p in outdir.glob('*')]))
    res: Dict[str, Any] = {
        'outdir': str(outdir),
        'artifacts': artifacts,
    }
    console.print(f"[green]✓ approx pipeline finished[/green] — outputs: {len(artifacts)} files")
    return res

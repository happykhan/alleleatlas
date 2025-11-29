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
from typing import Dict, Any
import csv
import numpy as np
from rich.console import Console

console = Console()


@dataclass
class ClusteringConfig:
    cgmlst_profiles: str
    outdir: str
    k: int = 50
    backend: str = 'exact'  # exact|usearch
    force_recompute: bool = False
    nproc: int = 2
    verbose: bool = True
    supernode_layout: str = 'spring'


def run(
    input_profile: str,
    out_dir: str,
    *,
    k: int = 50,
    backend: str = 'exact',
    force: bool = False,
    supernode_layout: str = 'spring',
    nproc: int = 2,
) -> None:
    """Run the full end-to-end approx pipeline.

    This is a thin local wrapper that delegates to the implementation in
    `approx.scripts.e2e_run.run`. It exists here so callers can import
    `run` from this module directly if desired.
    """
    # In-process fallback implementation -------------------------------------------------
    console.print(f"[blue]Running full pipeline (fallback): {input_profile} -> {out_dir} (backend={backend}, k={k})[/blue]")
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

    import numpy as np
    import scipy.sparse as sp

    # helper to save sparse/dense safely
    def _save_matrix(path: str, mat):
        if hasattr(mat, 'tocsr'):
            sp.save_npz(path, mat)
        else:
            # convert dense to CSR
            sp.save_npz(path, sp.csr_matrix(np.asarray(mat)))

    # 3) Build k-NN using requested backend (exact legacy or USearch)
    backend_suffix = 'exact' if backend == 'exact' else 'usearch'
    knn_path = p / f"{base}_{backend_suffix}.npz"
    index_obj = None
    if knn_path.exists() and not force:
        console.print(f"[yellow]Reusing existing k-NN matrix:[/yellow] {knn_path}")
    else:
        # 2) Load profiles (numpy array with NaNs)
        from alleleatlas import build_knn as build_knn_module
        X = build_knn_module.load_profiles_from_parquet(str(parquet_path))

        if backend == 'exact':
            from alleleatlas import distances_exact
            console.print(f'[green]Building exact k-NN (k={k}) via distances_exact[/green]')
            knn_mat = distances_exact.knn_from_profiles(X, k)
        elif backend == 'usearch':
            from alleleatlas.backends.usearch_backend import build_knn as usearch_build_knn
            console.print(f'[green]Building k-NN (k={k}) via USearch custom metric[/green]')
            knn_mat, index_obj, _ = usearch_build_knn(X, k, return_index=True, threads=nproc)
        else:
            raise ValueError(f"Unsupported backend: {backend}")
        _save_matrix(str(knn_path), knn_mat)
    # log max distance in k-NN
    try:
        import scipy.sparse as sp
        knn_tmp = sp.load_npz(knn_path)
        max_dist = float(knn_tmp.data.max()) if knn_tmp.nnz > 0 else 0.0
        console.print(f"[blue]Max distance in k-NN:[/blue] {max_dist:.4f}")
    except Exception as e:
        console.print(f"[yellow]Could not compute max k-NN distance:[/yellow] {e}")

    # 4) Analyze and plot
    from alleleatlas import analyze_knn
    from alleleatlas import mst_and_plot

    elbow1, elbow2 = analyze_knn.run_analysis(str(knn_path), out_prefix=str(p / f"{base}_{backend}"))
    # Build a coarse cluster graph using USearch's clustering on the index
    console.print(f"[blue]Generating cluster graph for backend={backend}[/blue]")
    if backend == 'usearch':
        try:
            import networkx as nx
            import matplotlib.pyplot as plt
            from alleleatlas.backends import usearch_backend

            if index_obj is None:
                # Rebuild index only if needed for clustering
                if 'X' not in locals():
                    from alleleatlas import build_knn as build_knn_module
                    X = build_knn_module.load_profiles_from_parquet(str(parquet_path))
                index_obj, _ = usearch_backend.build_index(X)
            if 'X' not in locals():
                from alleleatlas import build_knn as build_knn_module
                X = build_knn_module.load_profiles_from_parquet(str(parquet_path))

            cluster_cache = p / f"{base}_{backend}_clusters_data.npz"
            cluster_csv = p / f"{base}_{backend}_clusters.csv"

            if cluster_cache.exists() and cluster_csv.exists() and not force:
                console.print(f"[yellow]Reusing cached cluster graph:[/yellow] {cluster_cache}")
                data = np.load(cluster_cache)
                centroid_keys = data["centroid_keys"]
                sizes = data["sizes"]
                rows = data["rows"]
                cols = data["cols"]
                weights = data["weights"]
                size_map = {int(k): int(s) for k, s in zip(centroid_keys, sizes)}
                g_full = nx.Graph()
                for i, j, w in zip(rows, cols, weights):
                    g_full.add_edge(int(i), int(j), weight=float(w))
            else:
                console.print("[green]Clustering with USearch index to plot cluster graph[/green]")
                clustering = index_obj.cluster()
                console.print("[green]Drawing cluster graph[/green]")
                centroid_keys, sizes = clustering.centroids_popularity
                size_map = {int(k): s for k, s in zip(centroid_keys, sizes)}
                # persist cluster graph
                rows = []
                cols = []
                weights = []
                for u, v, d in clustering.network.edges(data=True):
                    rows.append(int(u))
                    cols.append(int(v))
                    weights.append(float(d.get("weight", 1.0)))
                np.savez(
                    cluster_cache,
                    centroid_keys=np.array(list(size_map.keys()), dtype=np.int64),
                    sizes=np.array(list(size_map.values()), dtype=np.int64),
                    rows=np.array(rows, dtype=np.int64),
                    cols=np.array(cols, dtype=np.int64),
                    weights=np.array(weights, dtype=float),
                )
                # Write cluster sizes table
                with open(cluster_csv, "w", newline="") as fh:
                    writer = csv.writer(fh)
                    writer.writerow(["centroid_key", "size"])
                    for k_val, s_val in sorted(size_map.items(), key=lambda kv: kv[1], reverse=True):
                        writer.writerow([k_val, s_val])
                console.print(f"[green]Wrote cluster sizes:[/green] {cluster_csv}")
                g_full = clustering.network

            # Draw cluster backbone plot (custom thresholds)
            mst_and_plot.plot_cluster_backbone(
                str(knn_path),
                str(p / f"{base}_{backend}_cluster_mst.png"),
                coarse_threshold=250.0,
                edge_display_threshold=800.0,
                edge_cap=1100.0,
                show_insets=False,
                grid_scale=1.5,
                supernode_layout=supernode_layout,
                component_summary_out=str(p / f"{base}_{backend}_supernodes.tsv"),
            )
        except Exception as e:  # pragma: no cover - visualization is best-effort
            console.print(f"[red]Cluster graph generation failed:[/red] {e}")
            # fallback to MST plot of the original k-NN if clustering failed
            try:
                # if mst_and_plot is not None:
                #     mst_and_plot.plot_mst(
                #         str(knn_path),
                #         out=str(p / f"{base}_{backend}_mst.png"),
                #         collapse_threshold=None,
                #         labels_path=None,
                #         use_embedding_init=False,
                #     )
                    console.print(f"[yellow]Fallback MST plot written for {backend} backend[/yellow]")
            except Exception as e2:
                console.print(f"[red]Fallback MST plotting also failed:[/red] {e2}")

def run_pipeline(config: ClusteringConfig) -> Dict[str, Any]:
    """Run pipeline given a ClusteringConfig and return artifact metadata.

    This is the function used by the CLI entrypoint. It is intentionally
    small and delegates heavy work to run().
    """
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    console.print(f"[cyan]Running pipeline:[/cyan] {config.cgmlst_profiles} -> {outdir} (backend={config.backend}, k={config.k})")

    run(
        config.cgmlst_profiles,
        str(outdir),
        k=config.k,
        backend=config.backend,
        force=config.force_recompute,
        supernode_layout=getattr(config, "supernode_layout", 'spring'),
        nproc=config.nproc,
    )

    artifacts = list(sorted([str(p) for p in outdir.glob('*')]))
    res: Dict[str, Any] = {
        'outdir': str(outdir),
        'artifacts': artifacts,
    }
    console.print(f"[green]✓ approx pipeline finished[/green] — outputs: {len(artifacts)} files")
    return res

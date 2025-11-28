"""Build k-NN graphs from normalized parquet profiles.

This module is a small dispatcher that loads a normalized parquet (as
produced by approx.convert_to_parquet) and delegates k-NN construction to a
selected backend implementation under ``approx.backends``.

See the backends directory for available implementations: hnsw, minhash_lsh,
and sklearn_nn. Each backend exposes a build_knn(profiles: np.ndarray, k: int)
function that returns a scipy sparse CSR matrix.

"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import typer
from rich.console import Console

console = Console()
app = typer.Typer()


def load_profiles_from_parquet(path: str) -> np.ndarray:
    """Load the parquet file and return an (n, d) numpy array of int64 profiles.

    This function expects the parquet produced by `convert_to_parquet` where
    ST is the first column and subsequent columns are integer loci. We return
    the loci-only matrix (drop ST column).
    """
    import pandas as pd
    df = pd.read_parquet(path)
    # assume first column is ST and drop it (it's an identifier, not a locus)
    if 'ST' in df.columns:
        loci = df.drop(columns=['ST'])
    else:
        loci = df.iloc[:, 1:]
    # convert to float to allow NaN for sentinel/missing values
    arr = loci.to_numpy(dtype=float)
    # detect sentinel large values and set to NaN (missing data)
    # sentinel heuristic: values much larger than typical allele values
    sentinel_thresh = 1e8
    arr[arr > sentinel_thresh] = np.nan
    return arr


@app.command()
def cli(
    profiles: str = typer.Argument(..., help='Input parquet with normalized profiles'),
    out: str = typer.Argument(..., help='Output .npz path for kNN (scipy.sparse npz)'),
    k: int = typer.Option(50, help='Number of neighbors'),
    backend: str = typer.Option('exact', help='Backend to use: exact|usearch|hnsw|minhash|sklearn'),
    usearch_dtype: str = typer.Option('f32', help='USearch dtype (f32/f16/f64/f8/i8/b1)'),
    # hnsw options
    hnsw_ef_construction: int = typer.Option(200, help='hnsw ef_construction'),
    hnsw_M: int = typer.Option(16, help='hnsw M parameter'),
    hnsw_ef_search: int = typer.Option(200, help='hnsw ef (search) parameter'),
    # minhash options
    num_perm: int = typer.Option(128, help='MinHash num_perm'),
    threshold: float = typer.Option(0.5, help='MinHash LSH threshold'),
):
    p = Path(out)
    p.parent.mkdir(parents=True, exist_ok=True)
    console.print(f"Loading profiles from [blue]{profiles}[/blue]")
    X = load_profiles_from_parquet(profiles)
    console.print(f"Profiles shape: {X.shape}")
    if backend == 'exact':
        from alleleatlas import distances_exact
    elif backend == 'usearch':
        from alleleatlas.backends.usearch_backend import build_knn as usearch_build_knn
    elif backend == 'hnsw':
        from alleleatlas.backends.hnsw import build_knn
    elif backend == 'minhash':
        from alleleatlas.backends.minhash_lsh import build_knn
    elif backend == 'sklearn':
        from alleleatlas.backends.sklearn_nn import build_knn
    else:
        console.print(f"[red]Unknown backend:[/red] {backend}")
        raise typer.Exit(code=1)

    console.print(f"Building k-NN (k={k}) using backend [green]{backend}[/green]")
    if backend == 'exact':
        knn_mat = distances_exact.knn_from_profiles(X, k)
    elif backend == 'usearch':
        knn_mat = usearch_build_knn(X, k, dtype=usearch_dtype)
    elif backend == 'hnsw':
        knn_mat = build_knn(X, k, ef_construction=hnsw_ef_construction, M=hnsw_M, ef_search=hnsw_ef_search)
    elif backend == 'minhash':
        knn_mat = build_knn(X, k, num_perm=num_perm, threshold=threshold)
    else:
        knn_mat = build_knn(X, k)
    import scipy.sparse as sp
    sp.save_npz(str(p), knn_mat)
    console.print(f"Wrote k-NN npz: {out}")


if __name__ == '__main__':
    app()

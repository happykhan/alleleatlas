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


def run(input_profile: str, out_dir: str, *, k: int = 50, backend: str = 'exact') -> None:
    """Run the full end-to-end approx pipeline.

    This is a thin local wrapper that delegates to the implementation in
    `approx.scripts.e2e_run.run`. It exists here so callers can import
    `run` from this module directly if desired.
    """
    # In-process fallback implementation -------------------------------------------------
    console.print(f"[blue]Running full pipeline (fallback): {input_profile} -> {out_dir} (backend={backend}, k={k})[/blue]")
    from pathlib import Path
    p = Path(out_dir)
    p.mkdir(parents=True, exist_ok=True)
    base = Path(input_profile).stem

    # 1) Convert input profile to parquet
    from alleleatlas import convert_to_parquet
    parquet = str(p / f"{base}.parquet")
    convert_to_parquet.convert_file(input_profile, parquet)

    # 2) Load profiles (numpy array with NaNs)
    from alleleatlas import build_knn as build_knn_module
    X = build_knn_module.load_profiles_from_parquet(parquet)

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
    if backend == 'exact':
        from alleleatlas import distances_exact
        console.print(f'[green]Building exact k-NN (k={k}) via distances_exact[/green]')
        knn_mat = distances_exact.knn_from_profiles(X, k)
        knn_path = str(p / f"{base}_exact.npz")
    elif backend == 'usearch':
        from alleleatlas.backends.usearch_backend import build_knn as usearch_build_knn
        console.print(f'[green]Building k-NN (k={k}) via USearch custom metric[/green]')
        knn_mat = usearch_build_knn(X, k)
        knn_path = str(p / f"{base}_usearch.npz")
    else:
        raise ValueError(f"Unsupported backend: {backend}")
    _save_matrix(knn_path, knn_mat)

    # 4) Analyze and plot
    from alleleatlas import analyze_knn
    from alleleatlas import mst_and_plot

    if analyze_knn is not None:
        analyze_knn.run_analysis(knn_path, out_prefix=str(p / f"{base}_{backend}"))

    # labels file heuristic: look for test_trees/{base}.tsv
    labels_path = None
    try:
        cand = Path('test_trees') / f"{base}.tsv"
        if cand.exists():
            labels_path = str(cand)
    except Exception:
        labels_path = None

    if mst_and_plot is not None:
        mst_and_plot.plot_mst(
            knn_path,
            out=str(p / f"{base}_{backend}_mst.png"),
            collapse_threshold=None,
            labels_path=labels_path,
            use_embedding_init=False,
        )

def run_pipeline(config: ClusteringConfig) -> Dict[str, Any]:
    """Run pipeline given a ClusteringConfig and return artifact metadata.

    This is the function used by the CLI entrypoint. It is intentionally
    small and delegates heavy work to run().
    """
    outdir = Path(config.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    console.print(f"[cyan]Running pipeline:[/cyan] {config.cgmlst_profiles} -> {outdir} (backend={config.backend}, k={config.k})")

    run(config.cgmlst_profiles, str(outdir), k=config.k, backend=config.backend)

    artifacts = list(sorted([str(p) for p in outdir.glob('*')]))
    res: Dict[str, Any] = {
        'outdir': str(outdir),
        'artifacts': artifacts,
    }
    console.print(f"[green]✓ approx pipeline finished[/green] — outputs: {len(artifacts)} files")
    return res

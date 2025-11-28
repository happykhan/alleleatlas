"""Small testbed / benchmark for alleleatlas pipelines.

This script is intentionally lightweight: it runs the exact distance
computation and (optionally) an approximate backend, times each step,
computes neighbor recall and MST edge overlap, and writes a JSON report.

Usage (from repo root):
    python -m tools.testbed --profile test_trees/paraC.profile --out out/testbed_run.json

It avoids heavy optional dependencies; if an approximate backend is not
available it will skip it and still report the exact timings.
"""
import argparse
import typer
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn
import json
import time
from pathlib import Path
import resource
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from alleleatlas import distances_exact
from alleleatlas.mst_exact import build_mst_from_distance_matrix
import numpy as np
from alleleatlas.backends import hamming, minhash_lsh, sklearn_nn, hnsw, rerank
from alleleatlas.backends import hybrid_union
import scipy.sparse as sp



def mem_kb():
    # peak RSS in KB (platform dependent but works on macOS/linux)
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_maxrss


def load_profiles_array(profile_path):
    # use convert_to_parquet.load_profiles if available; fallback to numpy load
    p = Path(profile_path)
    from alleleatlas.convert_to_parquet import load_profiles

    df, _ = load_profiles(str(p))
    # drop ST if present
    if 'ST' in df.columns:
        df = df.drop(columns=['ST'])
    arr = df.values.astype(float)
    return arr


def run_exact(profiles):
    t0 = time.time()
    legacy_arr = distances_exact.get_legacy_pairwise(profiles, method='numba')
    t1 = time.time()
    mem = mem_kb()
    D = legacy_arr[:, :, 0].astype(float)
    return {'time_s': t1 - t0, 'mem_kb': mem, 'shape': D.shape}


def run_exact_knn(profiles, k=10):
    t0 = time.time()
    knn = distances_exact.knn_from_profiles(profiles, k=k)
    t1 = time.time()
    return {'time_s': t1 - t0, 'knn': knn}



def recall_between_knn(reference_knn, approx_neighbors, k=10):
    # reference_knn: csr matrix
    ref = reference_knn.tocsr()
    n = ref.shape[0]
    total = 0
    hit = 0
    for i in range(n):
        ref_n = set(ref.indices[ref.indptr[i]: ref.indptr[i+1]])
        appr = set(approx_neighbors[i][:k])
        if len(ref_n) == 0:
            continue
        total += len(ref_n)
        hit += len(ref_n & appr)
    return hit / total if total > 0 else 0.0


def mst_edge_set_from_csr(knn_csr):
    # build full distance matrix by re-ranking candidates for MST
    # The MST builder expects a dense distance matrix; for speed we will
    # reconstruct a symmetric distance matrix where non-present edges are large.
    import numpy as np
    n = knn_csr.shape[0]
    INF = 10**9
    D = np.full((n, n), INF, dtype=float)
    for i in range(n):
        D[i, i] = 0.0
    rows = knn_csr.tocsr()
    for i in range(n):
        for idx in range(rows.indptr[i], rows.indptr[i+1]):
            j = rows.indices[idx]
            # symmetric placeholder; exact distance unknown here so set 1
            D[i, j] = 1.0
            D[j, i] = 1.0
    G = build_mst_from_distance_matrix(D)
    edges = set(tuple(sorted((u, v))) for u, v in G.edges())
    return edges


def neighbors_to_csr(neighbors, n):
    # neighbors: list[list[int]] per node
    rows = []
    cols = []
    data = []
    for i, nbrs in enumerate(neighbors):
        for j in nbrs:
            rows.append(i)
            cols.append(j)
            data.append(1)
    if len(rows) == 0:
        return sp.csr_matrix((n, n), dtype=int)
    mat = sp.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=int)
    return mat


def plot_summary(json_path, out_png):
    data = json.loads(Path(json_path).read_text(encoding='utf-8'))
    runs = data.get('runs', {})
    # filter to only backends that produced correctness metrics (exclude reference 'exact' rows)
    filtered = {k: v for k, v in runs.items() if ('recall' in v) or ('mst_edge_overlap' in v)}
    labels = list(filtered.keys())

    # collect metrics, allowing missing values
    recalls = [filtered[name].get('recall') if isinstance(filtered[name].get('recall'), (int, float)) else None for name in labels]
    overlaps = [filtered[name].get('mst_edge_overlap') if isinstance(filtered[name].get('mst_edge_overlap'), (int, float)) else None for name in labels]

    # plotting side-by-side bars for correctness metrics
    x = range(len(labels))
    width = 0.35
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4), sharey=False)

    # Recall plot
    ax = axes[0]
    recall_vals = [v if v is not None else 0.0 for v in recalls]
    missing_recall = [v is None for v in recalls]
    bars = ax.bar(x, recall_vals, width, color='C0')
    # dim bars where missing
    for i, m in enumerate(missing_recall):
        if m:
            bars[i].set_alpha(0.3)
    ax.set_ylim(0, 1)
    ax.set_title('Neighbor Recall')
    ax.set_ylabel('recall')

    # MST overlap plot
    ax2 = axes[1]
    overlap_vals = [v if v is not None else 0.0 for v in overlaps]
    missing_overlap = [v is None for v in overlaps]
    bars2 = ax2.bar(x, overlap_vals, width, color='C1')
    for i, m in enumerate(missing_overlap):
        if m:
            bars2[i].set_alpha(0.3)
    ax2.set_ylim(0, 1)
    ax2.set_title('MST edge overlap')
    ax2.set_ylabel('overlap')

    # rotate x labels under both plots
    for axx in axes:
        axx.set_xticks(x)
        axx.set_xticklabels(labels, rotation=90, ha='center')

    fig.subplots_adjust(bottom=0.35, wspace=0.4)
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches='tight')


app = typer.Typer()
console = Console()


@app.command()
def main(
    profile: str = typer.Option(..., help='Path to input profile (parquet/tsv)'),
    out: str = typer.Option(..., help='JSON output path'),
    k: int = typer.Option(10, help='k for k-NN')
):
    """Run the lightweight testbed/benchmark for alleleatlas backends."""
    console.print(f"[bold]Testbed:[/bold] loading profiles from [cyan]{profile}[/cyan]")
    profiles = load_profiles_array(profile)
    outd = {'input_rows': profiles.shape[0], 'input_cols': profiles.shape[1], 'runs': {}}

    with Progress(SpinnerColumn(), TextColumn('{task.description}'), transient=True) as prog:
        task = prog.add_task('Computing exact (legacy) pairwise distances...', start=False)
        prog.start_task(task)
        outd['runs']['exact'] = run_exact(profiles)
        prog.update(task, description='Computing exact k-NN...')
        ek = run_exact_knn(profiles, k=k)
        outd['runs']['exact_knn'] = {'time_s': ek['time_s']}
        exact_knn = ek['knn']

    console.print('[bold]Testbed:[/bold] benchmarking backends...')
    backends = [
        ('hamming', hamming),
        ('minhash_lsh', minhash_lsh),
        ('sklearn_nn', sklearn_nn),
        ('hnsw', hnsw),
        ('hybrid_union', hybrid_union),
        ('hnsw_rerank', rerank),
    ]

    for name, mod in backends:
        console.print(f' - Running backend [green]{name}[/green]')
        try:
            if hasattr(mod, 'build_knn'):
                t0 = time.time()
                knn = mod.build_knn(profiles, k)
                t1 = time.time()
                outd['runs'][name] = {'time_s': t1 - t0}
                if hasattr(knn, 'tocsr'):
                    recall = recall_between_knn(exact_knn, [list(knn.tocsr().indices[knn.tocsr().indptr[i]:knn.tocsr().indptr[i+1]]) for i in range(knn.shape[0])], k=k)
                    outd['runs'][name]['recall'] = recall
                    approx_csr = knn.tocsr()
                    n = approx_csr.shape[0]
                    rows, cols, data = [], [], []
                    for i in range(n):
                        cand = approx_csr.indices[approx_csr.indptr[i]: approx_csr.indptr[i+1]]
                        if cand.size == 0:
                            continue
                        vec = profiles[i]
                        dists = np.array([distances_exact._pair_distance_row(vec, profiles)[j] for j in cand], dtype=float)
                        for j, d in zip(cand, dists):
                            rows.append(i)
                            cols.append(int(j))
                            data.append(float(d))
                    reranked_csr = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
                    try:
                        exact_edges = mst_edge_set_from_csr(exact_knn)
                        approx_edges = mst_edge_set_from_csr(reranked_csr)
                        overlap = len(exact_edges & approx_edges) / len(exact_edges) if len(exact_edges) > 0 else 0.0
                    except Exception:
                        overlap = None
                    outd['runs'][name]['mst_edge_overlap'] = overlap
            elif name == 'hnsw_rerank' and hasattr(mod, 'hnsw_rerank'):
                t0 = time.time()
                knn = mod.hnsw_rerank(profiles, k)
                t1 = time.time()
                outd['runs'][name] = {'time_s': t1 - t0}
                recall = recall_between_knn(exact_knn, [list(knn.tocsr().indices[knn.tocsr().indptr[i]:knn.tocsr().indptr[i+1]]) for i in range(knn.shape[0])], k=k)
                outd['runs'][name]['recall'] = recall
                approx_edges = mst_edge_set_from_csr(knn)
                exact_edges = mst_edge_set_from_csr(exact_knn)
                overlap = len(exact_edges & approx_edges) / len(exact_edges) if len(exact_edges) > 0 else 0.0
                outd['runs'][name]['mst_edge_overlap'] = overlap
        except Exception as exc:
            console.print(f'[red]Skipped {name}:[/red] {exc}')
            outd['runs'][name] = {'skipped': True, 'reason': str(exc)}

    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(outd, indent=2), encoding='utf-8')
    console.print(f'[bold]Wrote[/bold] {out_path}')

    png_path = str(out_path.with_suffix('.png'))
    plot_summary(str(out_path), png_path)
    console.print(f'[bold]Wrote[/bold] {png_path}')


def _argparse_entrypoint():
    # keep backwards compatibility for scripts calling via argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--profile', required=True)
    parser.add_argument('--out', required=True, help='JSON output path')
    parser.add_argument('--k', type=int, default=10, help='k for k-NN')
    parser.add_argument('--minhash-perm', type=int, default=128, help='MinHash permutations')
    parser.add_argument('--use-legacy', action='store_true', help='Use legacy numba getDistance implementation as ground truth')
    args = parser.parse_args()
    # delegate to typer-style function directly (only pass supported args)
    main(profile=args.profile, out=args.out, k=args.k)


if __name__ == '__main__':
    app()

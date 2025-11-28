"""Benchmark k-NN backend implementations on synthetic data.

Produces simple timing and memory usage reports for each backend implementation.
"""
from __future__ import annotations

import time
import numpy as np
from rich.console import Console

try:
    import psutil
except Exception:
    psutil = None

console = Console()


def gen_profiles(n: int, d: int, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.integers(0, 1000, size=(n, d), dtype=np.int64)


def run_backend(backend_name: str, X: np.ndarray, k: int, **kwargs):
    console.print(f"Running backend [green]{backend_name}[/green] n={X.shape[0]} d={X.shape[1]} k={k}")
    if backend_name == 'hnsw':
        from approx.backends.hnsw import build_knn
    elif backend_name == 'minhash':
        from approx.backends.minhash_lsh import build_knn
    elif backend_name == 'sklearn':
        from approx.backends.sklearn_nn import build_knn
    else:
        raise RuntimeError('unknown backend')
    # memory snapshot before
    p = psutil.Process() if psutil else None
    mem_before = p.memory_info().rss if p else None
    t0 = time.perf_counter()
    M = build_knn(X, k, **kwargs)
    t1 = time.perf_counter()
    mem_after = p.memory_info().rss if p else None
    mem_delta = (mem_after - mem_before) if (mem_before and mem_after) else None
    console.print(f"Time: {t1 - t0:.3f}s, nnz={M.nnz}")
    if mem_delta is not None:
        console.print(f"Memory delta (RSS): {mem_delta / 1024**2:.2f} MiB")


def main():
    # small experiments
    sizes = [(100, 64), (500, 64)]
    backends = ['hnsw', 'minhash', 'sklearn']
    for n, d in sizes:
        X = gen_profiles(n, d)
        for b in backends:
            try:
                if b == 'hnsw':
                    run_backend(b, X, k=10, ef_construction=200, M=16, ef_search=200)
                elif b == 'minhash':
                    run_backend(b, X, k=10, num_perm=128, threshold=0.5)
                else:
                    run_backend(b, X, k=10)
            except Exception as e:
                console.print(f"[red]Backend {b} failed:[/red] {e}")


if __name__ == '__main__':
    main()
"""Benchmark k-NN backend implementations on synthetic data.

Produces simple timing and memory usage reports for each backend implementation.
"""
from __future__ import annotations

import time
import numpy as np
from rich.console import Console

console = Console()


def gen_profiles(n: int, d: int, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.integers(0, 1000, size=(n, d), dtype=np.int64)


def run_backend(backend_name: str, X: np.ndarray, k: int):
    console.print(f"Running backend [green]{backend_name}[/green] n={X.shape[0]} d={X.shape[1]} k={k}")
    if backend_name == 'hnsw':
        from approx.backends.hnsw import build_knn
    elif backend_name == 'minhash':
        from approx.backends.minhash_lsh import build_knn
    elif backend_name == 'sklearn':
        from approx.backends.sklearn_nn import build_knn
    else:
        raise RuntimeError('unknown backend')
    t0 = time.perf_counter()
    M = build_knn(X, k)
    t1 = time.perf_counter()
    console.print(f"Time: {t1 - t0:.3f}s, nnz={M.nnz}")


def main():
    # small experiments
    sizes = [(100, 64), (500, 64)]
    backends = ['hnsw', 'minhash', 'sklearn']
    for n, d in sizes:
        X = gen_profiles(n, d)
        for b in backends:
            try:
                run_backend(b, X, k=10)
            except Exception as e:
                console.print(f"[red]Backend {b} failed:[/red] {e}")


if __name__ == '__main__':
    main()

"""USearch backend for cgMLST distances.

This backend uses USearch with a custom metric that matches the
missing-aware allelic distance (count of differing loci where both
profiles have non-missing values; missing coded as 0). It can build a
USearch index and derive a k-NN sparse matrix compatible with the rest
of the pipeline.

Note: we keep the metric as a Python callable (Metric) rather than a
numba-compiled function to avoid crashes seen with compiled pointers on
some macOS/Python combos. This trades a small amount of speed for
stability.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def _to_int_profiles(profiles: np.ndarray) -> np.ndarray:
    """Convert float profiles with NaN-missing to int32 with 0-missing."""
    arr = profiles.astype(float)
    arr[np.isnan(arr)] = 0.0
    return arr.astype(np.int32, copy=False)


def _build_metric():
    """Return a USearch Metric for cgMLST allelic distance."""
    try:
        from usearch.index import MetricKind, MetricSignature, CompiledMetric
    except ImportError as e:  # pragma: no cover - optional dep
        raise ImportError("USearch python bindings are required for backend='usearch'") from e

    def cgmlst_distance(a, b):
        diff = 0
        comparable = 0
        for ai, bi in zip(a, b):
            if ai == 0 or bi == 0:
                continue
            comparable += 1
            if ai != bi:
                diff += 1
        if comparable == 0:
            return float(len(a))
        return float(diff)

    return CompiledMetric(
        cgmlst_distance,
        kind=MetricKind.L2sq,  # nominal; real distance is provided
        signature=MetricSignature.ArrayArray,
    )


def build_index(profiles: np.ndarray, dtype: str = "f32"):
    """Build and return a USearch Index along with the keys used."""
    try:
        from usearch.index import Index
    except ImportError as e:  # pragma: no cover - optional dep
        raise ImportError("USearch python bindings are required for backend='usearch'") from e

    mat = _to_int_profiles(profiles)
    # Restrict to f32 to align with the compiled metric signature
    if dtype != 'f32':
        raise ValueError("USearch backend currently supports dtype='f32' only for the custom metric")
    mat = np.ascontiguousarray(mat.astype(np.float32))

    n, d = mat.shape
    metric = _build_metric()
    index = Index(ndim=d, metric=metric, dtype=dtype)
    keys = np.arange(n, dtype=np.int64)
    index.add(keys, mat)
    return index, keys


def build_knn(profiles: np.ndarray, k: int, dtype: str = "f32") -> sp.csr_matrix:
    """Build k-NN sparse matrix using USearch index search."""
    index, keys = build_index(profiles, dtype=dtype)
    n = len(keys)
    if n == 0 or k <= 0:
        return sp.csr_matrix((n, n), dtype=float)
    rows = []
    cols = []
    data = []
    mat = _to_int_profiles(profiles)
    k_eff = min(k, n - 1)
    for i in range(n):
        query = mat[i]
        res = index.search(query, k_eff + 1)  # +1 to allow self, filter out below
        for key, dist in zip(res.keys, res.distances):
            if key == keys[i]:
                continue
            rows.append(i)
            cols.append(int(key))
            data.append(float(dist))
    if not rows:
        return sp.csr_matrix((n, n), dtype=float)
    return sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()

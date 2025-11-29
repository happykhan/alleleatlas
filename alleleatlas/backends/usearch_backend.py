"""USearch backend for cgMLST distances.

This backend uses USearch with a custom metric that matches the
missing-aware allelic distance (count of differing loci where both
profiles have non-missing values; missing coded as 0). It can build a
USearch index and derive a k-NN sparse matrix compatible with the rest
of the pipeline.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from functools import lru_cache


def _to_int_profiles(profiles: np.ndarray) -> np.ndarray:
    """Convert float profiles with NaN-missing to int32 with 0-missing."""
    arr = profiles.astype(float)
    arr[np.isnan(arr)] = 0.0
    return arr.astype(np.int32, copy=False)


@lru_cache(maxsize=None)
def _build_metric(ndim: int, allowed_missing: float = 0.05):
    """Return a USearch Metric for cgMLST allelic distance matching legacy behavior."""
    try:
        from usearch.index import MetricKind, MetricSignature, CompiledMetric
        from numba import cfunc, types, carray
    except ImportError as e:  # pragma: no cover - optional dep
        raise ImportError("USearch python bindings are required for backend='usearch'") from e

    signature = types.float32(types.CPointer(types.float32), types.CPointer(types.float32))
    am = types.float32(allowed_missing)
    n_loci = types.int32(ndim)

    @cfunc(signature)
    def cgmlst_distance(a_ptr, b_ptr):
        a = carray(a_ptr, ndim)
        b = carray(b_ptr, ndim)
        ql = types.float32(0.0)
        for i in range(ndim):
            if a[i] > 0:
                ql += 1.0

        rl = types.float32(0.0)
        ad = types.float32(1e-4)
        al = types.float32(1e-4)
        for i in range(ndim):
            bi = b[i]
            if bi > 0:
                rl += 1.0
                ai = a[i]
                if ai > 0:
                    al += 1.0
                    if ai != bi:
                        ad += 1.0

        ll = max(ql, rl) - am * n_loci
        ll2 = ql - am * n_loci

        if ll2 > al:
            ad += ll2 - al
            al = ll2

        # legacy v1 (not returned but kept for parity)
        _v1 = ad / al * n_loci + 0.5  # noqa: F841

        if ll > al:
            ad += ll - al
            al = ll
        v0 = ad / al * n_loci + 0.5
        return types.float32(types.int32(v0))

    return CompiledMetric(
        pointer=cgmlst_distance.address,
        kind=MetricKind.L2sq,  # nominal; real distance is provided
        signature=MetricSignature.ArrayArray,
    )


def build_index(profiles: np.ndarray, dtype: str = "f32", allowed_missing: float = 0.05):
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
    metric = _build_metric(d, float(allowed_missing))
    index = Index(ndim=d, metric=metric, dtype=dtype)
    keys = np.arange(n, dtype=np.int64)
    index.add(keys, mat)
    return index, keys


def build_knn(
    profiles: np.ndarray,
    k: int,
    dtype: str = "f32",
    return_index: bool = False,
    allowed_missing: float = 0.05,
    exact: bool = False,
    threads: int = 0,
):
    """Build k-NN sparse matrix using USearch index search.

    If `return_index` is True, also return the underlying USearch index and keys.
    """
    index, keys = build_index(profiles, dtype=dtype, allowed_missing=allowed_missing)
    n = len(keys)
    if n == 0 or k <= 0:
        empty = sp.csr_matrix((n, n), dtype=float)
        return (empty, index, keys) if return_index else empty
    rows = []
    cols = []
    data = []
    mat = np.ascontiguousarray(_to_int_profiles(profiles).astype(np.float32))
    k_eff = min(k, n - 1)
    for i in range(n):
        query = mat[i]
        res = index.search(query, k_eff + 1, exact=exact, threads=threads)  # +1 to allow self, filter out below
        for key, dist in zip(res.keys, res.distances):
            if key == keys[i]:
                continue
            rows.append(i)
            cols.append(int(key))
            data.append(float(dist))
    if not rows:
        knn = sp.csr_matrix((n, n), dtype=float)
    else:
        knn = sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    return (knn, index, keys) if return_index else knn

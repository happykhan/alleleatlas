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


from functools import lru_cache
import numpy as np

@lru_cache(maxsize=None)
def _build_metric_hclink(
    ndim: int,
    allowed_missing: float = 0.03,   # use 0.03 for HierCC parity
    max_gaps: int | None = None,
):
    """Return a USearch Metric implementing a HIERCC-style allelic distance.

    Distance:
      - counts differing loci where both profiles are non-missing
      - handles missing data with a query_core/common_core trade-off
      - caps distance at n_loci when too many shared gaps or extreme mismatch
    """
    try:
        from usearch.index import MetricKind, MetricSignature, CompiledMetric
        from numba import cfunc, types, carray
    except ImportError as e:
        raise ImportError("USearch python bindings are required for backend='usearch'") from e

    # USearch Index(ArrayArray) will always pass ndim=ndim, so we capture it in the closure.
    n_loci = ndim
    am = np.float32(allowed_missing)
    # If max_gaps is None, disable that cap by setting it larger than possible.
    max_gaps_val = np.int32(max_gaps if max_gaps is not None else n_loci + 1)
    n_f32 = np.float32(n_loci)

    signature = types.float32(
        types.CPointer(types.float32),
        types.CPointer(types.float32),
    )

    @cfunc(signature)
    def hiercc_distance(a_ptr, b_ptr):
        a = carray(a_ptr, n_loci)
        b = carray(b_ptr, n_loci)

        gaps_a = 0
        gaps_b = 0
        gaps_both = 0
        distance = 0

        for i in range(n_loci):
            av = a[i]
            bv = b[i]
            if av == bv:
                if av == 0:
                    gaps_both += 1
            else:
                if av != 0 and bv != 0:
                    distance += 1
                if av == 0:
                    gaps_a += 1
                if bv == 0:
                    gaps_b += 1

        # Cap distance if too many shared gaps or distance already maximal
        if gaps_both >= max_gaps_val or distance >= n_loci:
            return np.float32(n_loci)

        # Perfect match with no gaps in either profile
        if distance == 0 and gaps_a == 0 and gaps_b == 0:
            return np.float32(0.0)

        # HIERCC-style distance
        query_core  = n_f32 - np.float32(gaps_a + gaps_both) - am * n_f32
        common_core = n_f32 - np.float32(gaps_a + gaps_b + gaps_both)

        if common_core >= query_core:
            if common_core == 0.0:
                return n_f32
            cc_distance = (n_f32 * np.float32(distance)) / common_core + np.float32(0.5)
        else:
            if query_core == 0.0:
                return n_f32
            cc_distance = (
                n_f32 * np.float32(distance + (query_core - common_core)) / query_core
                + np.float32(0.5)
            )

        if cc_distance > n_f32:
            cc_distance = n_f32
        # Return as float32, but integral-valued in practice
        return np.float32(cc_distance)

    return CompiledMetric(
        pointer=hiercc_distance.address,
        kind=MetricKind.Unknown,           # not a true L2; Unknown is more honest
        signature=MetricSignature.ArrayArray,
    )



@lru_cache(maxsize=None)
def _build_metric_hiercc(ndim: int, allowed_missing: float = 0.03):
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


def build_index(
    profiles: np.ndarray,
    dtype: str = "f32",
    allowed_missing: float = 0.05,
    connectivity: int = 16,       # NEW: expose connectivity
    expansion_add: int = 128,     # NEW
    expansion_search: int = 256,  # NEW
    metric_mode: str = "hclink",  # NEW: choose metric mode
):
    """Build and return a USearch Index along with the keys used."""
    try:
        from usearch.index import Index
    except ImportError as e:  # pragma: no cover - optional dep
        raise ImportError("USearch python bindings are required for backend='usearch'") from e

    mat = _to_int_profiles(profiles)
    # Restrict to f32 to align with the compiled metric signature
    if dtype != "f32":
        raise ValueError("USearch backend currently supports dtype='f32' only for the custom metric")
    mat = np.ascontiguousarray(mat.astype(np.float32))
    n, d = mat.shape
    if metric_mode == "hclink":
        metric = _build_metric_hclink(d, float(allowed_missing))
    else:
        metric = _build_metric_hiercc(d, float(allowed_missing))
    # NEW: configure index more aggressively for higher recall
    index = Index(
        ndim=d,
        metric=metric,
        dtype=dtype,
        connectivity=connectivity,
        expansion_add=expansion_add,
        expansion_search=expansion_search,
    )
    # (setting expansion_search again is harmless but explicit)
    index.expansion_search = expansion_search  # NEW

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
    connectivity: int = 16,       # NEW: plumb through to build_index
    expansion_add: int = 128,     # NEW
    expansion_search: int = 256,  # NEW
    oversample_factor: float = 2.0,  # NEW: ask for more neighbours than we keep
    # advanced option:
    rerank_exact: bool = False,   # NEW: small exact rerank on candidates
    metric_mode: str = "hclink",  # NEW: choose metric mode
):
    """Build k-NN sparse matrix using USearch index search.

    If `return_index` is True, also return the underlying USearch index and keys.
    """
    index, keys = build_index(
        profiles,
        dtype=dtype,
        allowed_missing=allowed_missing,
        connectivity=connectivity,
        expansion_add=expansion_add,
        expansion_search=expansion_search,
    )
    n = len(keys)
    if n == 0 or k <= 0:
        empty = sp.csr_matrix((n, n), dtype=float)
        return (empty, index, keys) if return_index else empty

    mat = np.ascontiguousarray(_to_int_profiles(profiles).astype(np.float32))
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []

    k_eff = min(k, n - 1)
    # NEW: search for more neighbours than we keep, then truncate
    if oversample_factor < 1.0:
        oversample_factor = 1.0
    k_search = min(int(np.ceil(k_eff * oversample_factor)), n - 1)

    # optional: metric reuse for reranking
    d = mat.shape[1]
    if metric_mode == "hclink":
        metric = _build_metric_hclink(d, float(allowed_missing))
    else:
        metric = _build_metric_hiercc(d, float(allowed_missing))
    try:
        from usearch.index import Index as UIndex  # for rerank
    except ImportError:
        rerank_exact = False  # silently disable if not available

    for i in range(n):
        query = mat[i]

        # +1 to allow self, filter out below
        res = index.search(query, k_search + 1, exact=exact, threads=threads)

        # Filter out self and keep up to k_search candidates
        cand_keys = []
        cand_dists = []
        for key, dist in zip(res.keys, res.distances):
            if key == keys[i]:
                continue
            cand_keys.append(int(key))
            cand_dists.append(float(dist))
            if len(cand_keys) >= k_search:
                break

        # Optional: build a tiny exact rerank index on candidates
        if rerank_exact and not exact and cand_keys:
            cand_keys_arr = np.asarray(cand_keys, dtype=np.int64)
            cand_profiles = mat[cand_keys_arr]

            rerank_index = UIndex(
                ndim=d,
                metric=metric,
                dtype=dtype,
                connectivity=0,
                expansion_add=0,
                expansion_search=0,
            )
            rerank_index.add(cand_keys_arr, cand_profiles, copy=False)

            reranked = rerank_index.search(query, k_eff, exact=True)
            use_keys = [int(k_) for k_ in reranked.keys]
            use_dists = [float(d_) for d_ in reranked.distances]
        else:
            # Just truncate approximate results to k_eff
            use_keys = cand_keys[:k_eff]
            use_dists = cand_dists[:k_eff]

        for key_j, dist_j in zip(use_keys, use_dists):
            rows.append(i)
            cols.append(key_j)
            data.append(dist_j)

    if not rows:
        knn = sp.csr_matrix((n, n), dtype=float)
    else:
        knn = sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()

    return (knn, index, keys) if return_index else knn

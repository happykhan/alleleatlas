"""Utilities to rerank candidate neighbors (HNSW candidate generation + exact rerank).

This module provides a small helper that builds an HNSW index to obtain a
candidate set for each query and then re-ranks those candidates using exact
Euclidean distances on the original numeric vectors. The function ensures the
query index is never returned as its own neighbor.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from typing import Optional


def hnsw_rerank(
    profiles: np.ndarray,
    k: int,
    *,
    candidate_mult: int = 6,
    ef_construction: int = 500,
    M: int = 32,
    ef_search: int = 500,
    verbose: bool = False,
):
    """Generate candidates with HNSW then rerank by exact L2.

    Parameters
    - profiles: (n, d) numeric array
    - k: final neighbor count
    - candidate_mult: multiply k to produce candidate pool (k_cand = k * candidate_mult)
    - ef_construction, M, ef_search: hnswlib parameters
    - verbose: print progress

    Returns a scipy.sparse.csr_matrix (n, n) with distances for top-k neighbors.
    """
    try:
        import hnswlib
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError("hnswlib not installed") from exc

    n, d = profiles.shape
    X = profiles.astype(np.float32)

    idx = hnswlib.Index(space='l2', dim=d)
    idx.init_index(max_elements=n, ef_construction=ef_construction, M=M)
    idx.add_items(X, np.arange(n))
    idx.set_ef(ef_search)

    k_cand = min(n, int(k * candidate_mult))
    if verbose:
        print(f'hnsw_rerank: querying k_cand={k_cand} (n={n}, d={d})')
    labels, _ = idx.knn_query(X, k=k_cand)

    rows = []
    cols = []
    data = []

    # exact rerank per query, ensuring we exclude the query index itself
    for i in range(n):
        cand = np.unique(labels[i])
        # exclude query itself if present
        cand = cand[cand != i]
        if cand.size == 0:
            continue
        vec = profiles[i].astype(np.float64)
        cvecs = profiles[cand].astype(np.float64)
        dists = np.linalg.norm(cvecs - vec, axis=1)
        if dists.size > k:
            order = np.argpartition(dists, k)[:k]
            order = order[np.argsort(dists[order])]
            sel = cand[order]
            seld = dists[order]
        else:
            order = np.argsort(dists)
            sel = cand[order][:k]
            seld = dists[order][:k]

        rows.extend([i] * len(sel))
        cols.extend(sel.tolist())
        data.extend(seld.tolist())

    mat = sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    if verbose:
        print(f'hnsw_rerank: built csr nnz={mat.nnz}')
    return mat

"""Hybrid union backend: MinHash + Faiss binary/HNSW candidates, then unioned.

This module provides a convenience function `build_knn` that collects candidates
from multiple backends and returns a CSR matrix of candidates (not reranked).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from typing import Optional

from alleleatlas.backends import minhash_lsh
try:
    from alleleatlas.backends import faiss_binary
except Exception:
    faiss_binary = None
try:
    from alleleatlas.backends import hnsw as hnsw_backend
except Exception:
    hnsw_backend = None


def build_knn(profiles: np.ndarray, k: int, *, minhash_perms: int = 256, minhash_k: int = 50,
              faiss_dim_bits: int = 16384, use_faiss_ivf: bool = False, faiss_nprobe: int = 4,
              hnsw_params: Optional[dict] = None) -> sp.csr_matrix:
    """Collect candidate neighbors from MinHash and Faiss/HNSW and union them.

    Returns CSR where entries indicate candidate membership (value 1.0).
    """
    n = profiles.shape[0]
    candidate_sets = [set() for _ in range(n)]

    # MinHash candidates
    try:
        mh = minhash_lsh.build_knn(profiles, minhash_k, perms=minhash_perms)
        rows = mh.tocsr()
        for i in range(n):
            candidate_sets[i].update(rows.indices[rows.indptr[i]: rows.indptr[i+1]])
    except Exception:
        # ignore minhash failures
        pass

    # Faiss binary candidates
    if faiss_binary is not None:
        try:
            fb = faiss_binary.build_knn(profiles, k=min(minhash_k, n-1), dim_bits=faiss_dim_bits, use_ivf=use_faiss_ivf, nprobe=faiss_nprobe)
            rows = fb.tocsr()
            for i in range(n):
                candidate_sets[i].update(rows.indices[rows.indptr[i]: rows.indptr[i+1]])
        except Exception:
            pass

    # HNSW fallback
    if hnsw_backend is not None and hnsw_params is not None:
        try:
            hb = hnsw_backend.build_knn(profiles if profiles.ndim==2 else profiles.astype(float), k=minhash_k, **hnsw_params)
            rows = hb.tocsr() if hasattr(hb, 'tocsr') else hb
            for i in range(n):
                candidate_sets[i].update(rows.indices[rows.indptr[i]: rows.indptr[i+1]])
        except Exception:
            pass

    # Build CSR
    rows = []
    cols = []
    data = []
    for i in range(n):
        for j in candidate_sets[i]:
            if j == i:
                continue
            rows.append(i)
            cols.append(int(j))
            data.append(1.0)
    if len(rows) == 0:
        return sp.csr_matrix((n, n), dtype=float)
    return sp.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=float)

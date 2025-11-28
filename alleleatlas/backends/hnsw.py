"""HNSW backend using hnswlib to build approximate k-NN for integer profiles.

This module exposes ``build_knn(profiles, k, *, ef_construction, M, ef_search, verbose)``
that returns a :class:`scipy.sparse.csr_matrix` containing pairwise neighbor
distances for each row.

Notes:
- Inputs are integer-valued profiles (n, d). We convert to float for the
  hnswlib index. For large d you may want to pre-project or normalize first.
"""
from __future__ import annotations

import numpy as np


def build_knn(
    profiles: np.ndarray,
    k: int,
    *,
    ef_construction: int = 200,
    M: int = 16,
    ef_search: int = 200,
    verbose: bool = False,
):
    """Build approximate k-NN using hnswlib.

    Parameters
    - profiles: integer-valued array (n, d)
    - k: number of neighbors
    - ef_construction, M: hnswlib index construction parameters
    - ef_search: ef value set on the index for querying
    - verbose: if True, print progress via Rich-friendly prints
    """
    try:
        import hnswlib
    except Exception as exc:  # pragma: no cover - optional deps
        raise RuntimeError('hnswlib not installed') from exc

    n, d = profiles.shape
    # convert to float vectors for hnsw (simple cast). Users may prefer
    # normalization, PCA, or MinHash pre-processing for high dimensional data.
    X = profiles.astype(float)

    idx = hnswlib.Index(space='l2', dim=d)
    idx.init_index(max_elements=n, ef_construction=ef_construction, M=M)
    idx.add_items(X, np.arange(n))
    idx.set_ef(ef_search)
    labels, distances = idx.knn_query(X, k=k + 1)
    # drop self-match (first column)
    labels = labels[:, 1:]
    distances = distances[:, 1:]

    import scipy.sparse as sp
    rows = np.repeat(np.arange(n), k)
    cols = labels.reshape(-1)
    data = distances.reshape(-1)
    Mtx = sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    if verbose:
        print(f'hnsw: built approx kNN n={n} d={d} k={k} nnz={Mtx.nnz}')
    return Mtx

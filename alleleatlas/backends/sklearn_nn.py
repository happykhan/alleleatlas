"""Simple backend using scikit-learn NearestNeighbors for exact/approx NN.

This implementation is intended as a reliable baseline for small-to-moderate
datasets and for testing/benchmarking.
"""
from __future__ import annotations

import numpy as np


def build_knn(profiles: np.ndarray, k: int):
    try:
        from sklearn.neighbors import NearestNeighbors
    except Exception as exc:  # pragma: no cover - optional dependency
        raise RuntimeError('scikit-learn not available') from exc
    n, d = profiles.shape
    if n == 0:
        import scipy.sparse as sp
        return sp.csr_matrix((0, 0))
    X = profiles.astype(float)
    nn = NearestNeighbors(n_neighbors=min(k+1, n), algorithm='auto', metric='euclidean')
    nn.fit(X)
    distances, labels = nn.kneighbors(X, n_neighbors=min(k+1, n))
    # drop self-match if present at position 0
    if labels.shape[1] > 0 and (labels[:, 0] == np.arange(n)).all():
        labels = labels[:, 1:]
        distances = distances[:, 1:]
    import scipy.sparse as sp
    k_eff = labels.shape[1]
    rows = np.repeat(np.arange(n), k_eff)
    cols = labels.reshape(-1)
    data = distances.reshape(-1)
    M = sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    return M

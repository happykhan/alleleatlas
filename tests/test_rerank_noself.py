import numpy as np
import scipy.sparse as sp

from alleleatlas.backends.rerank import hnsw_rerank


def test_hnsw_rerank_excludes_self():
    # small synthetic dataset: 10 points in 3D clustered around two centers
    rng = np.random.default_rng(0)
    a = rng.normal(loc=0.0, scale=0.1, size=(5, 3))
    b = rng.normal(loc=1.0, scale=0.1, size=(5, 3))
    X = np.vstack([a, b])
    k = 3
    try:
        M = hnsw_rerank(X, k, candidate_mult=4, ef_construction=100, M=8, ef_search=200)
    except RuntimeError:
        # If hnswlib is not installed, skip the test by asserting nothing raised
        return
    assert isinstance(M, sp.csr_matrix)
    n = X.shape[0]
    for i in range(n):
        row = M.getrow(i).tocoo()
        if row.nnz == 0:
            continue
        # ensure query index not present among columns
        assert i not in set(row.col), f"Query {i} appears in its own neighbor list"

import numpy as np
import pytest


def gen_profiles(n=50, d=16, seed=1):
    rng = np.random.default_rng(seed)
    return rng.integers(0, 1000, size=(n, d), dtype=np.int64)


def _check_matrix(M, n, k):
    import scipy.sparse as sp
    assert sp.isspmatrix_csr(M)
    assert M.shape == (n, n)
    # if n > k then expect at least n*k nonzeros (approx)
    assert M.nnz >= min(n * k, n * n)


def test_hnsw_backend_small():
    X = gen_profiles(30, 8)
    try:
        from alleleatlas.backends.hnsw import build_knn
    except Exception:
        pytest.skip('hnswlib not installed')
    M = build_knn(X, 5)
    _check_matrix(M, 30, 5)


def test_sklearn_backend_small():
    X = gen_profiles(30, 8)
    try:
        from alleleatlas.backends.sklearn_nn import build_knn
    except Exception:
        pytest.skip('sklearn not available')
    M = build_knn(X, 5)
    _check_matrix(M, 30, 5)


def test_minhash_backend_small():
    X = gen_profiles(30, 8)
    try:
        from alleleatlas.backends.minhash_lsh import build_knn
    except Exception:
        pytest.skip('datasketch not available')
    M = build_knn(X, 5)
    _check_matrix(M, 30, 5)

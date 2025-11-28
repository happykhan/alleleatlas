import numpy as np
from alleleatlas.distances_exact import (
    pairwise_legacy_numba,
    pairwise_legacy_python,
)


def compare_arrays(a, b):
    # both are int arrays shaped (n,n,2)
    assert a.shape == b.shape
    assert np.array_equal(a[:, :, 0], b[:, :, 0]), "mismatch counts differ"
    assert np.array_equal(a[:, :, 1], b[:, :, 1]), "comparable counts differ"


def test_parity_random_small_seeds():
    rng = np.random.default_rng(42)
    seeds = list(range(5))
    for s in seeds:
        rng = np.random.default_rng(s)
        N = 12
        D = 40
        profiles = rng.integers(1, 20, size=(N, D)).astype(float)
        a = pairwise_legacy_numba(profiles, allowed_missing=0.05)
        b = pairwise_legacy_python(profiles, allowed_missing=0.05)
        compare_arrays(a, b)


def test_parity_various_shapes():
    rng = np.random.default_rng(7)
    shapes = [(4, 8), (8, 16), (16, 32)]
    for N, D in shapes:
        profiles = rng.integers(1, 50, size=(N, D)).astype(float)
        a = pairwise_legacy_numba(profiles, allowed_missing=0.0)
        b = pairwise_legacy_python(profiles, allowed_missing=0.0)
        compare_arrays(a, b)


def test_parity_edge_cases():
    # all identical rows
    profiles = np.ones((6, 10), dtype=float)
    a = pairwise_legacy_numba(profiles)
    b = pairwise_legacy_python(profiles)
    compare_arrays(a, b)

    # sequential rows
    profiles = np.vstack([np.arange(1, 11) + i for i in range(6)]).astype(float)
    a = pairwise_legacy_numba(profiles)
    b = pairwise_legacy_python(profiles)
    compare_arrays(a, b)

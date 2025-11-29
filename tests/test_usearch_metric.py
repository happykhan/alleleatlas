import numpy as np
import pytest


def test_usearch_metric_matches_legacy():
    pytest.importorskip("usearch")
    from alleleatlas.backends import usearch_backend
    from alleleatlas.distances_exact import pairwise_legacy_numba

    profiles = np.array(
        [
            [1, 2, 0, 3],  # missing at idx 2
            [1, 2, 4, 3],  # one mismatch at idx 2 (missing vs 4)
            [1, 5, 4, 3],  # two mismatches (idx1, idx2)
            [0, 0, 0, 0],  # all missing
        ],
        dtype=float,
    )

    legacy = pairwise_legacy_numba(profiles)[:, :, 0].astype(float)

    # build index and exact distances via search
    index, keys = usearch_backend.build_index(profiles)
    usearch_dists = np.zeros_like(legacy)
    mat = np.ascontiguousarray(profiles)
    for i, key in enumerate(keys):
        res = index.search(mat[i], len(keys), exact=True)
        for k, dist in zip(res.keys, res.distances):
            usearch_dists[i, int(k)] = dist

    # self distances should be zero in legacy; set same for comparison
    np.fill_diagonal(usearch_dists, 0)

    assert np.allclose(legacy, usearch_dists)

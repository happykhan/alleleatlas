import numpy as np
from alleleatlas import distances_exact


def test_pairwise_small_triangle():
    # three profiles with single locus each
    # A:0, B:1, C:1 -> distances: AB=1, AC=1, BC=0
    profiles = np.array([[0], [1], [1]], dtype=float)
    D = distances_exact.pairwise_distances(profiles)
    # symmetric, zero diagonal
    assert D.shape == (3, 3)
    assert D[0, 0] == 0
    assert D[1, 1] == 0
    assert D[2, 2] == 0
    assert D[0, 1] == 1
    assert D[0, 2] == 1
    assert D[1, 2] == 0


def test_knn_from_profiles_basic():
    # 4 profiles, two clusters: [0,0],[0,0],[1,1],[1,1]
    profiles = np.array([[0, 0], [0, 0], [1, 1], [1, 1]], dtype=float)
    # expect k=1 neighbours to be the other item in the same cluster
    k = 1
    knn = distances_exact.knn_from_profiles(profiles, k=k)
    # knn is scipy csr matrix with shape (4, 4)
    assert knn.shape == (4, 4)
    # for each row, check that the single neighbor is in the expected index set
    rows = knn.tocsr()
    for i in range(4):
        neighbors = rows.indices[rows.indptr[i]: rows.indptr[i+1]]
        assert len(neighbors) == k
        nb = neighbors[0]
        if i in (0, 1):
            assert nb in (0, 1) and nb != i
        else:
            assert nb in (2, 3) and nb != i

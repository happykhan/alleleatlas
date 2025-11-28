import numpy as np
from alleleatlas.mst_exact import build_mst_from_distance_matrix, mst_to_newick


def test_mst_simple_triangle():
    D = np.array([[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.0]])
    G = build_mst_from_distance_matrix(D)
    assert G.number_of_edges() == 2
    nwk = mst_to_newick(G)
    assert nwk.endswith(';')

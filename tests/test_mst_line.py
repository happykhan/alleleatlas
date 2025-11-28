import numpy as np
import networkx as nx
from scipy.spatial import distance
from alleleatlas.mst_and_plot import build_component_mst


def test_mst_on_line_points():
    # deterministic 1D points on a line (distinct coordinates)
    xs = np.linspace(0, 10, 12)
    X = xs.reshape(-1, 1).astype(float)
    # exact dense distances
    D = distance.cdist(X, X, metric='euclidean')
    # set diagonal to zero for sparse representation
    np.fill_diagonal(D, 0.0)
    # reference MST: chain connecting consecutive points
    G_ref = nx.from_numpy_array(D)
    T_ref = nx.minimum_spanning_tree(G_ref, weight='weight')

    # pipeline: use exact dense distances as the k-NN adjacency matrix (dense -> sparse)
    import scipy.sparse as sp
    n = X.shape[0]
    knn = sp.csr_matrix(D)
    # convert k-NN sparse to MST using the same helper
    G_pipeline = build_component_mst(knn)
    T_pipe = nx.minimum_spanning_tree(G_pipeline)

    # compare edge sets (unordered)
    def edge_set(T):
        return set(frozenset((u, v)) for u, v in T.edges())

    Eref = edge_set(T_ref)
    Epipe = edge_set(T_pipe)

    # For points on a line the MST should be exactly the chain (n-1 edges)
    assert len(Eref) == n - 1
    assert len(Epipe) == n - 1
    # exact equality of edge sets
    assert Eref == Epipe

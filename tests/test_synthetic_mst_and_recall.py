import math
import random
import numpy as np
import networkx as nx
import scipy.spatial as sps
import scipy.sparse as sp
from alleleatlas.mst_and_plot import build_component_mst
from alleleatlas.backends.sklearn_nn import build_knn as sklearn_build
from alleleatlas.backends.rerank import hnsw_rerank


def exact_knn(X, k):
    D = sps.distance.cdist(X, X, metric='euclidean')
    np.fill_diagonal(D, np.inf)
    idx = np.argpartition(D, k, axis=1)[:, :k]
    return [list(row) for row in idx], D


def recall_at_k(ref_idx, pred_idx, k):
    n = len(ref_idx)
    tot = 0.0
    for i in range(n):
        r = set(ref_idx[i][:k])
        p = set(pred_idx[i][:k])
        if len(r) == 0:
            continue
        tot += len(r & p) / len(r)
    return tot / n


def edge_set(T):
    return set(frozenset((u, v)) for u, v in T.edges())


def mst_from_distance_matrix(D):
    G = nx.from_numpy_array(D)
    return nx.minimum_spanning_tree(G, weight='weight')


def test_ring_points_mst():
    # points sampled uniformly on a circle; MST is ring minus the longest chord
    n = 20
    theta = np.linspace(0, 2 * math.pi, n, endpoint=False)
    X = np.column_stack([np.cos(theta), np.sin(theta)])
    D = sps.distance.cdist(X, X, metric='euclidean')
    np.fill_diagonal(D, 0.0)

    T_ref = mst_from_distance_matrix(D)
    # pipeline: build k-NN as dense adjacency and compute MST
    knn = sp.csr_matrix(D)
    G_pipeline = build_component_mst(knn)
    T_pipe = nx.minimum_spanning_tree(G_pipeline)

    # check MST edges count and that exactly one of the longest ring chords is missing
    assert len(edge_set(T_ref)) == n - 1
    assert len(edge_set(T_pipe)) == n - 1
    # edge overlap should be near 1.0 for exact distances
    overlap = len(edge_set(T_ref) & edge_set(T_pipe)) / (n - 1)
    assert overlap == 1.0


def test_grid_points_mst():
    # 2D grid points; MST should connect nearest neighbors consistently
    xs = np.arange(4)
    ys = np.arange(3)
    xv, yv = np.meshgrid(xs, ys)
    X = np.column_stack([xv.ravel(), yv.ravel()])
    D = sps.distance.cdist(X, X, metric='euclidean')
    np.fill_diagonal(D, 0.0)

    T_ref = mst_from_distance_matrix(D)
    knn = sp.csr_matrix(D)
    G_pipeline = build_component_mst(knn)
    T_pipe = nx.minimum_spanning_tree(G_pipeline)

    # Expect near identical MST for exact distances
    overlap = len(edge_set(T_ref) & edge_set(T_pipe)) / (len(edge_set(T_ref)) or 1)
    assert overlap == 1.0


def test_two_blobs_single_bridge():
    # two well-separated Gaussian blobs; MST should have exactly one inter-cluster edge
    rng = np.random.default_rng(1)
    a = rng.normal(loc=0.0, scale=0.1, size=(10, 2))
    b = rng.normal(loc=5.0, scale=0.1, size=(10, 2))
    X = np.vstack([a, b])
    D = sps.distance.cdist(X, X, metric='euclidean')
    np.fill_diagonal(D, 0.0)

    labels = np.array([0] * 10 + [1] * 10)

    T_ref = mst_from_distance_matrix(D)
    # count inter-cluster edges in MST
    inter = 0
    for u, v in T_ref.edges():
        if labels[u] != labels[v]:
            inter += 1
    # ideally exactly one bridge between clusters
    assert inter == 1


def _random_tree(n_tips, seed=0):
    # generate random bifurcating tree with random branch lengths
    rng = random.Random(seed)
    # create a list of current nodes (tips or internal)
    nodes = list(range(n_tips))
    next_id = n_tips
    edges = []  # (u, v, length)
    while len(nodes) > 1:
        i = rng.randrange(len(nodes))
        j = rng.randrange(len(nodes))
        if i == j:
            continue
        u = nodes.pop(max(i, j))
        v = nodes.pop(min(i, j))
        # create new parent
        parent = next_id
        next_id += 1
        # assign branch lengths
        lu = rng.random() * 0.5 + 0.01
        lv = rng.random() * 0.5 + 0.01
        edges.append((parent, u, lu))
        edges.append((parent, v, lv))
        nodes.append(parent)
    # return edge list
    return edges


def _patristic_matrix(n_tips, edges):
    # build adjacency and compute pairwise path lengths between tips
    # nodes are numbered 0..(n_tips-1) for tips and higher for internals
    ids = set()
    for u, v, length in edges:
        ids.add(u)
        ids.add(v)
    G = nx.Graph()
    for u, v, length in edges:
        G.add_edge(u, v, weight=length)
    # compute distances between tips
    D = np.zeros((n_tips, n_tips), dtype=float)
    for i in range(n_tips):
        lengths = nx.single_source_dijkstra_path_length(G, i, weight='weight')
        for j in range(n_tips):
            D[i, j] = lengths.get(j, 0.0)
    return D


def test_simulated_phylogeny_mst_recovery():
    # simulate a small phylogeny, evolve binary loci, and check MST overlap
    n_tips = 12
    edges = _random_tree(n_tips, seed=2)
    D_true = _patristic_matrix(n_tips, edges)
    # reference MST on true patristic distances
    T_ref = mst_from_distance_matrix(D_true)

    # simulate binary loci: probability of flip per branch ~ 1 - exp(-mu * length)
    L = 500
    mu = 1.0
    rng = np.random.default_rng(3)
    # build tree graph
    G = nx.Graph()
    for u, v, length in edges:
        G.add_edge(u, v, weight=length)

    # root at last internal node (max node id)
    root = max(u for u, v, _ in edges)
    # perform post-order traversal and simulate
    states = {root: np.zeros(L, dtype=np.int8)}
    # build parent-child map by orienting edges away from root
    Tnx = nx.bfs_tree(G, root)
    for parent in nx.topological_generations(Tnx):
        for p in parent:
            for child in Tnx.successors(p):
                bl = G[p][child]['weight']
                p_flip = 1.0 - math.exp(-mu * bl)
                flips = rng.random(L) < p_flip
                states[child] = np.bitwise_xor(states[p], flips.astype(np.int8))

    # extract tip vectors
    tip_vecs = np.vstack([states[i] for i in range(n_tips)])
    # Hamming distances
    H = sps.distance.cdist(tip_vecs, tip_vecs, metric='hamming') * L
    np.fill_diagonal(H, 0.0)

    T_sim = mst_from_distance_matrix(H)

    overlap = len(edge_set(T_ref) & edge_set(T_sim)) / (n_tips - 1)
    # we expect moderate recovery (not perfect), require >50% edges recovered
    assert overlap >= 0.5


def test_recall_between_exact_and_backends():
    # small Gaussian dataset
    rng = np.random.default_rng(10)
    X = rng.normal(size=(40, 5))
    k = 5
    ref_idx, D = exact_knn(X, k)

    # sklearn backend
    M = sklearn_build((X * 100).astype(int), k)
    # convert sklearn sparse to neighbor lists (sorted by distance)
    neigh = []
    for i in range(M.shape[0]):
        row = M.getrow(i).tocoo()
        if row.nnz == 0:
            neigh.append([])
            continue
        idxs = row.col[np.argsort(row.data)]
        neigh.append(list(idxs))
    r_sk = recall_at_k(ref_idx, neigh, k)
    assert r_sk >= 0.99

    # optional: reranked hnsw if available
    try:
        rer = hnsw_rerank(X.astype(float), k, candidate_mult=6)
    except Exception:
        return
    neigh_r = []
    for i in range(rer.shape[0]):
        row = rer.getrow(i).tocoo()
        if row.nnz == 0:
            neigh_r.append([])
            continue
        idxs = row.col[np.argsort(row.data)]
        neigh_r.append(list(idxs))
    r_h = recall_at_k(ref_idx, neigh_r, k)
    assert r_h >= 0.95

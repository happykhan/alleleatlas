#!/usr/bin/env python3
"""Consume a .knn.npz, cluster the graph, compute medoids and MST between clusters.

Outputs:
 - <out_prefix>.clusters.tsv  (sample_id \t cluster_id)
 - <out_prefix>.medoids.tsv   (cluster_id \t sample_id)
 - <out_prefix>.cluster_mst.tsv (u_cluster \t v_cluster \t weight)
 - <out_prefix>.summary.json  (basic stats)

Usage:
  python3 scripts/pipeline_cluster_and_mst.py --knn out.knn.npz --names names.txt --out out_prefix

"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import minimum_spanning_tree

try:
    import igraph as ig
    import leidenalg
    HAVE_LEIDEN = True
except Exception:
    HAVE_LEIDEN = False

try:
    import community as community_louvain
    HAVE_LOUVAIN = True
except Exception:
    HAVE_LOUVAIN = False


def read_names(names_path: Path):
    with open(names_path, 'r') as fh:
        return [l.strip() for l in fh]


def run_leiden_on_csr(csr: sp.csr_matrix, resolution: float = 1.0):
    # build igraph from csr
    sources, targets = csr.nonzero()
    weights = csr.data
    g = ig.Graph(n=csr.shape[0], edges=list(zip(sources.tolist(), targets.tolist())), directed=False)
    g.es['weight'] = weights.tolist()
    partition = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, weights='weight', resolution_parameter=resolution)
    labels = np.array([p for p in partition.membership], dtype=int)
    return labels


def run_louvain_on_csr(csr: sp.csr_matrix):
    G = csr.tocsr()
    # networkx conversion for louvain is expensive; build minimal adjacency dict
    import networkx as nx
    Gnx = nx.from_scipy_sparse_matrix(G, edge_attribute='weight')
    part = community_louvain.best_partition(Gnx, weight='weight')
    labels = np.array([part[i] for i in range(G.shape[0])], dtype=int)
    return labels


def compute_medoid(cluster_ids: np.ndarray, dist_mem: np.memmap, n:int):
    # dist_mem is condensed distances; for medoid choose sample with minimal sum distance to others in cluster
    medoids = {}
    unique = np.unique(cluster_ids)
    for cid in unique:
        idxs = np.where(cluster_ids == cid)[0]
        if len(idxs) == 1:
            medoids[int(cid)] = int(idxs[0])
            continue
        # compute pairwise distances within idxs by reading condensed memmap
        m = len(idxs)
        D = np.zeros((m, m), dtype=float)
        for a in range(m):
            i = idxs[a]
            for b in range(a):
                j = idxs[b]
                # condensed index: ensure i>j order for function
                if i > j:
                    ii, jj = j, i
                else:
                    ii, jj = i, j
                # compute condensed idx
                idx = n * ii - (ii * (ii + 1) // 2) + (jj - ii - 1)
                d = float(dist_mem[idx])
                D[a, b] = d
                D[b, a] = d
        sums = D.sum(axis=1)
        medoids[int(cid)] = int(idxs[np.argmin(sums)])
    return medoids


def cluster_and_mst(knn_path: Path, names_path: Path, out_prefix: Path, use_leiden: bool = True):
    knn = sp.load_npz(str(knn_path))
    n = knn.shape[0]
    names = read_names(names_path)
    assert len(names) == n

    if HAVE_LEIDEN and use_leiden:
        labels = run_leiden_on_csr(knn)
    elif HAVE_LOUVAIN:
        labels = run_louvain_on_csr(knn)
    else:
        # Fallback: use connected components on the symmetric k-NN graph
        print('Warning: No Leiden/Louvain available; falling back to connected components')
        from scipy.sparse.csgraph import connected_components
        sym = knn.tocsr()
        sym = (sym + sym.transpose()) / 2.0
        n_components, labels = connected_components(sym, directed=False, return_labels=True)

    # write cluster assignments
    out_clusters = out_prefix.with_suffix('.clusters.tsv')
    with open(out_clusters, 'w') as fh:
        for i, lab in enumerate(labels):
            fh.write(f"{names[i]}\t{lab}\n")

    # compute medoids
    # try to find dist memmap from sibling path
    dist0_guess = str(knn_path).replace('.knn.npz', '.dist0.memmap')
    dist_mem = np.memmap(dist0_guess, dtype=np.int32, mode='r')
    medoids = compute_medoid(labels, dist_mem, n)
    out_medoids = out_prefix.with_suffix('.medoids.tsv')
    with open(out_medoids, 'w') as fh:
        for cid, midx in sorted(medoids.items()):
            fh.write(f"{cid}\t{names[midx]}\n")

    # compute cluster-level MST: create cluster centroid distance matrix (by medoid distances)
    clusters = sorted(medoids.keys())
    C = len(clusters)
    if C <= 1:
        print('Only one cluster found; skipping cluster MST')
        return
    # build full CxC distance matrix using medoid pairwise distances
    cent = np.zeros((C, C), dtype=float)
    for a in range(C):
        for b in range(a+1, C):
            ia = medoids[clusters[a]]
            ib = medoids[clusters[b]]
            if ia > ib:
                ii, jj = ib, ia
            else:
                ii, jj = ia, ib
            idx = n * ii - (ii * (ii + 1) // 2) + (jj - ii - 1)
            cent[a, b] = float(dist_mem[idx])
            cent[b, a] = cent[a, b]

    # MST on centroids
    # convert to sparse
    coo = sp.coo_matrix(cent)
    mst = minimum_spanning_tree(coo)
    mst_coo = mst.tocoo()
    out_mst = out_prefix.with_suffix('.cluster_mst.tsv')
    with open(out_mst, 'w') as fh:
        for u, v, w in zip(mst_coo.row, mst_coo.col, mst_coo.data):
            fh.write(f"{clusters[u]}\t{clusters[v]}\t{float(w)}\n")

    summary = {
        'n_samples': n,
        'n_clusters': C,
        'clusters': C,
    }
    open(str(out_prefix.with_suffix('.summary.json')), 'w').write(json.dumps(summary))
    print('Wrote cluster outputs:', out_clusters, out_medoids, out_mst)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--knn', required=True, help='Path to .knn.npz')
    p.add_argument('--names', required=True, help='Path to names.txt produced by condensed step')
    p.add_argument('--out', required=True, help='Output prefix')
    args = p.parse_args()
    cluster_and_mst(Path(args.knn), Path(args.names), Path(args.out))


if __name__ == '__main__':
    main()

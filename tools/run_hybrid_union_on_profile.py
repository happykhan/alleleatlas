#!/usr/bin/env python3
"""Run hybrid union backend (MinHash + Faiss/HNSW) on profile and rerank candidates.

Writes JSON to out/hybrid_union_verify.json
"""
from pathlib import Path
import json
import time
from alleleatlas.backends import hybrid_union
import numpy as np
import scipy.sparse as sp
from alleleatlas import distances_exact
from alleleatlas.mst_exact import build_mst_from_distance_matrix

OUT = Path('out')
OUT.mkdir(exist_ok=True)


def load_profiles_array(profile_path):
    from alleleatlas.convert_to_parquet import load_profiles
    df, _ = load_profiles(profile_path)
    if 'ST' in df.columns:
        df = df.drop(columns=['ST'])
    return df.values.astype(int)


def mst_edge_set_from_csr(knn_csr):
    n = knn_csr.shape[0]
    INF = 10**9
    D = np.full((n, n), INF, dtype=float)
    for i in range(n):
        D[i, i] = 0.0
    rows = knn_csr.tocsr()
    for i in range(n):
        for idx in range(rows.indptr[i], rows.indptr[i+1]):
            j = rows.indices[idx]
            D[i, j] = 1.0
            D[j, i] = 1.0
    G = build_mst_from_distance_matrix(D)
    edges = set(tuple(sorted((u, v))) for u, v in G.edges())
    return edges


def rerank_candidates_with_exact(profiles, approx_csr, K=5):
    n = approx_csr.shape[0]
    rows, cols, data = [], [], []
    for i in range(n):
        cand = approx_csr.indices[approx_csr.indptr[i]: approx_csr.indptr[i+1]]
        if cand.size == 0:
            continue
        vec = profiles[i]
        dists = np.array([distances_exact._pair_distance_row(vec, profiles)[j] for j in cand], dtype=float)
        if len(dists) > K:
            sel = np.argpartition(dists, K)[:K]
            chosen = [cand[s] for s in sel]
            chosen_d = dists[sel]
        else:
            chosen = list(cand)
            chosen_d = dists
        for j, d in zip(chosen, chosen_d):
            rows.append(i)
            cols.append(int(j))
            data.append(float(d))
    if len(rows) == 0:
        return sp.csr_matrix((n, n), dtype=float)
    return sp.csr_matrix((data, (rows, cols)), shape=(n, n))


def run(profile='test_trees/paraC.profile', out='out/hybrid_union_verify.json', k=5):
    profiles = load_profiles_array(profile)
    n = profiles.shape[0]

    # exact baseline
    legacy = distances_exact.get_legacy_pairwise(profiles, method='numba')
    exact_D = legacy[:, :, 0].astype(float)
    exact_knn = np.argsort(exact_D, axis=1)[:, 1:k+1]
    exact_csr = sp.csr_matrix(([1]*(n*k), (np.repeat(np.arange(n), k), exact_knn.reshape(-1))), shape=(n,n))
    exact_edges = mst_edge_set_from_csr(exact_csr)

    cfg = {'minhash_perms': 256, 'minhash_k': 50, 'faiss_dim_bits': 16384, 'use_faiss_ivf': False}
    t0 = time.time()
    approx = hybrid_union.build_knn(profiles, k=k, minhash_perms=cfg['minhash_perms'], minhash_k=cfg['minhash_k'], faiss_dim_bits=cfg['faiss_dim_bits'], use_faiss_ivf=cfg['use_faiss_ivf'])
    t1 = time.time()

    rows = approx.tocsr()
    # compute recall_direct (count hit if any approx candidate contains the reference neighbor)
    total = 0
    hit = 0
    for i in range(n):
        ref = set(exact_knn[i])
        appr = set(rows.indices[rows.indptr[i]: rows.indptr[i+1]])
        total += len(ref)
        hit += len(ref & appr)
    recall_direct = hit/total if total>0 else 0.0

    # rerank to final k
    reranked = rerank_candidates_with_exact(profiles, approx, K=k)
    approx_edges = mst_edge_set_from_csr(reranked)
    overlap = len(exact_edges & approx_edges)/len(exact_edges) if len(exact_edges)>0 else 0.0

    results = {'input_rows': n, 'runs': [{'cfg': cfg, 'status': 'done', 'build_time': t1-t0, 'recall_direct': recall_direct, 'recall_reranked': None, 'mst_edge_overlap': overlap}]}
    Path(out).write_text(json.dumps(results, indent=2))
    print('Wrote', out)


if __name__ == '__main__':
    run()

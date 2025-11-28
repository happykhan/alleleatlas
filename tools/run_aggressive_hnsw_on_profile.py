#!/usr/bin/env python3
"""Run an aggressive one-hot HNSW config (large D, M, ef) and rerank.

Writes JSON to out/hnsw_onehot_aggressive.json
"""
from pathlib import Path
import json
import time
import hashlib
import numpy as np
import scipy.sparse as sp
from alleleatlas.backends import hnsw
from alleleatlas import distances_exact
from alleleatlas.mst_exact import build_mst_from_distance_matrix

OUT = Path('out')
OUT.mkdir(exist_ok=True)


def one_hot_hash_encode(profiles, dim=65536, hash_name='blake2b'):
    n, d = profiles.shape
    X = np.zeros((n, dim), dtype=float)
    hfn = hashlib.blake2b if hash_name=='blake2b' else hashlib.md5
    for i in range(n):
        for j in range(d):
            v = profiles[i, j]
            if isinstance(v, float) and np.isnan(v):
                continue
            key = f"{j}_{int(v)}"
            h = hfn(key.encode('utf8')).hexdigest()
            idx = int(h, 16) % dim
            X[i, idx] = 1.0
    return X


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


def run(profile='test_trees/paraC.profile', out='out/hnsw_onehot_aggressive.json', k=5):
    profiles = load_profiles_array(profile)
    n = profiles.shape[0]

    legacy = distances_exact.get_legacy_pairwise(profiles, method='numba')
    exact_D = legacy[:, :, 0].astype(float)
    exact_knn = np.argsort(exact_D, axis=1)[:, 1:k+1]
    exact_csr = sp.csr_matrix(([1]*(n*k), (np.repeat(np.arange(n), k), exact_knn.reshape(-1))), shape=(n,n))
    exact_edges = mst_edge_set_from_csr(exact_csr)

    cfg = {'dim': 65536, 'M': 64, 'ef_construction': 500, 'ef_search': 2000, 'candidate_k': 50}
    Xh = one_hot_hash_encode(profiles, dim=cfg['dim'], hash_name='blake2b')
    t0 = time.time()
    knn_h = hnsw.build_knn(Xh, min(n-1, cfg['candidate_k']), ef_construction=cfg['ef_construction'], M=cfg['M'], ef_search=cfg['ef_search'])
    t1 = time.time()
    rows = knn_h.tocsr()

    total = 0
    hit = 0
    for i in range(n):
        ref = set(exact_knn[i])
        appr = set(rows.indices[rows.indptr[i]: rows.indptr[i+1]])
        total += len(ref)
        hit += len(ref & appr)
    recall_direct = hit/total if total>0 else 0.0

    reranked = rerank_candidates_with_exact(profiles, knn_h, K=k)
    approx_edges = mst_edge_set_from_csr(reranked)
    overlap = len(exact_edges & approx_edges)/len(exact_edges) if len(exact_edges)>0 else 0.0

    results = {'input_rows': n, 'runs': [{'cfg': cfg, 'status': 'done', 'build_time': t1-t0, 'recall_direct': recall_direct, 'recall_reranked': None, 'mst_edge_overlap': overlap}]}
    Path(out).write_text(json.dumps(results, indent=2))
    print('Wrote', out)


if __name__ == '__main__':
    run()

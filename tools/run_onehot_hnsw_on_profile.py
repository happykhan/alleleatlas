#!/usr/bin/env python3
"""Run one-hot HNSW experiments on a real profile.
Writes JSON to out/hnsw_onehot_verify.json
"""
import json
import time
import hashlib
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from alleleatlas.backends import hnsw, minhash_lsh
from alleleatlas import distances_exact
from alleleatlas.mst_exact import build_mst_from_distance_matrix

OUT = Path('out')
OUT.mkdir(exist_ok=True)


def one_hot_hash_encode(profiles, dim=4096, hash_name='md5'):
    n, d = profiles.shape
    X = np.zeros((n, dim), dtype=float)
    for i in range(n):
        for j in range(d):
            v = profiles[i, j]
            if isinstance(v, float) and np.isnan(v):
                continue
            key = f"{j}_{int(v)}"
            if hash_name == 'md5':
                h = hashlib.md5(key.encode('utf8')).hexdigest()
            else:
                h = hashlib.blake2b(key.encode('utf8')).hexdigest()
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


def run(profile='test_trees/paraC.profile', out='out/hnsw_onehot_verify.json'):
    profiles = load_profiles_array(profile)
    n = profiles.shape[0]
    # exact baseline
    legacy = distances_exact.get_legacy_pairwise(profiles, method='numba')
    exact_D = legacy[:, :, 0].astype(float)
    k = 5
    exact_knn = np.argsort(exact_D, axis=1)[:, 1:k+1]
    exact_csr = sp.csr_matrix(([1]*(n*k), (np.repeat(np.arange(n), k), exact_knn.reshape(-1))), shape=(n,n))
    exact_edges = mst_edge_set_from_csr(exact_csr)

    # configs to test: top-by-recalltime plus a high-recall config
    topf = Path('out/hnsw_sweep_top_by_recalltime.json')
    top = json.loads(topf.read_text())['top_by_recalltime']
    extra = {'dim': 16384, 'M': 32, 'ef_construction': 200, 'ef_search': 200, 'candidate_multiplier': 3, 'candidate_k': 15}

    dims = [4096, 8192, 16384]
    results = {'input_rows': n, 'runs': []}

    for dim in dims:
        Xh = one_hot_hash_encode(profiles, dim=dim, hash_name='md5')
        for cfg in top + [extra]:
            cfg = dict(cfg)
            cfg['dim'] = dim
            k_cand = min(n-1, cfg.get('candidate_k', 5))
            # build HNSW on hashed encoding
            t0 = time.time()
            try:
                knn_h = hnsw.build_knn(Xh, k_cand, ef_construction=cfg['ef_construction'], M=cfg['M'], ef_search=cfg['ef_search'])
            except Exception as exc:
                results['runs'].append({'cfg': cfg, 'status': 'error', 'error': str(exc)})
                continue
            t1 = time.time()
            build_time = t1 - t0
            # extract neighbors
            if hasattr(knn_h, 'tocsr'):
                approx_neighbors = [list(knn_h.tocsr().indices[knn_h.tocsr().indptr[i]:knn_h.tocsr().indptr[i+1]]) for i in range(n)]
            else:
                approx_neighbors = knn_h
            # compute recall between exact top-k and approx candidates
            total = 0
            hit = 0
            for i in range(n):
                ref = set(exact_knn[i])
                appr = set(approx_neighbors[i][:k])
                total += len(ref)
                hit += len(ref & appr)
            recall_direct = hit/total if total>0 else 0.0

            # rerank to final k and compute MST overlap
            reranked = rerank_candidates_with_exact(profiles, knn_h, K=k)
            approx_edges = mst_edge_set_from_csr(reranked)
            overlap = len(exact_edges & approx_edges)/len(exact_edges) if len(exact_edges)>0 else 0.0

            results['runs'].append({
                'cfg': cfg,
                'status': 'done',
                'build_time': build_time,
                'nnz': int(knn_h.nnz) if hasattr(knn_h, 'nnz') else None,
                'recall_direct': recall_direct,
                'recall_reranked': None, # could compute by producing final knn lists after rerank
                'mst_edge_overlap': overlap,
            })
    Path(out).write_text(json.dumps(results, indent=2))
    print('Wrote', out)

if __name__ == '__main__':
    run()

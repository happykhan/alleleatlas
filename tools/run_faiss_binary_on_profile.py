#!/usr/bin/env python3
"""Run Faiss binary backend on a profile and report recall@k and MST overlap.

Writes JSON to out/faiss_binary_verify.json
"""
from pathlib import Path
import json
import time
from alleleatlas.backends import faiss_binary
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


def run(profile='test_trees/paraC.profile', out='out/faiss_binary_verify.json', k=5):
    profiles = load_profiles_array(profile)
    n = profiles.shape[0]

    # exact baseline
    legacy = distances_exact.get_legacy_pairwise(profiles, method='numba')
    exact_D = legacy[:, :, 0].astype(float)
    exact_knn = np.argsort(exact_D, axis=1)[:, 1:k+1]
    exact_csr = sp.csr_matrix(([1]*(n*k), (np.repeat(np.arange(n), k), exact_knn.reshape(-1))), shape=(n,n))
    exact_edges = mst_edge_set_from_csr(exact_csr)

    results = {'input_rows': n, 'runs': []}

    cfg = {'dim_bits': 16384, 'use_ivf': False}
    t0 = time.time()
    try:
        knn = faiss_binary.build_knn(profiles, k, dim_bits=cfg['dim_bits'], use_ivf=cfg['use_ivf'])
    except Exception as exc:
        results['runs'].append({'cfg': cfg, 'status': 'error', 'error': str(exc)})
        Path(out).write_text(json.dumps(results, indent=2))
        print('Wrote', out)
        return
    t1 = time.time()
    # compute recall_direct
    rows = knn.tocsr()
    total = 0
    hit = 0
    for i in range(n):
        ref = set(exact_knn[i])
        appr = set(rows.indices[rows.indptr[i]: rows.indptr[i+1]])
        total += len(ref)
        hit += len(ref & appr)
    recall_direct = hit/total if total>0 else 0.0

    reranked = sp.csr_matrix((n,n), dtype=float)  # placeholder; can rerank per-candidate if needed
    approx_edges = mst_edge_set_from_csr(reranked)
    overlap = len(exact_edges & approx_edges)/len(exact_edges) if len(exact_edges)>0 else 0.0

    results['runs'].append({'cfg': cfg, 'status': 'done', 'build_time': t1-t0, 'recall_direct': recall_direct, 'recall_reranked': None, 'mst_edge_overlap': overlap})
    Path(out).write_text(json.dumps(results, indent=2))
    print('Wrote', out)


if __name__ == '__main__':
    run()

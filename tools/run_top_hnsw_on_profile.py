#!/usr/bin/env python3
"""Run top HNSW configs (from out/hnsw_sweep_top_by_recalltime.json) on a real profile and measure recall/MST overlap.

Usage:
    .venv/bin/python tools/run_top_hnsw_on_profile.py --profile test_trees/paraC.profile --out out/hnsw_top_verify.json
"""
import json
from pathlib import Path
import time
import numpy as np
import scipy.sparse as sp
from alleleatlas import distances_exact
from alleleatlas.backends import hnsw
from alleleatlas.mst_exact import build_mst_from_distance_matrix


def load_profiles_array(profile_path):
    from alleleatlas.convert_to_parquet import load_profiles
    df, _ = load_profiles(profile_path)
    if 'ST' in df.columns:
        df = df.drop(columns=['ST'])
    return df.values.astype(int)


def total_time(entry):
    return entry.get('build_time',0)+entry.get('rerank_time',0)+entry.get('minhash_time',0)


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


def rerank_candidates_with_exact(profiles, approx_csr):
    n = approx_csr.shape[0]
    rows, cols, data = [], [], []
    for i in range(n):
        cand = approx_csr.indices[approx_csr.indptr[i]: approx_csr.indptr[i+1]]
        if cand.size == 0:
            continue
        vec = profiles[i]
        # compute exact distances for these candidates using distances_exact helper
        dists = np.array([distances_exact._pair_distance_row(vec, profiles)[j] for j in cand], dtype=float)
        for j, d in zip(cand, dists):
            rows.append(i); cols.append(int(j)); data.append(float(d))
    reranked = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    return reranked


def main(profile='test_trees/paraC.profile', out='out/hnsw_top_verify.json'):
    topf = Path('out/hnsw_sweep_top_by_recalltime.json')
    if not topf.exists():
        raise SystemExit('Missing out/hnsw_sweep_top_by_recalltime.json; run tools/analyze_hnsw_sweep.py first')
    top = json.loads(topf.read_text())['top_by_recalltime']
    profiles = load_profiles_array(profile)
    n = profiles.shape[0]
    results = {'input_rows': n, 'configs': []}
    # compute exact knn and edges
    legacy = distances_exact.get_legacy_pairwise(profiles, method='numba')
    exact_D = legacy[:, :, 0].astype(float)
    # naive exact knn from exact_D
    import numpy as np
    k = 5
    exact_knn = np.argsort(exact_D, axis=1)[:, 1:k+1]
    exact_csr = sp.csr_matrix(([1]* (n*k), (np.repeat(np.arange(n), k), exact_knn.reshape(-1))), shape=(n,n))
    exact_edges = mst_edge_set_from_csr(exact_csr)

    for ci, cfg in enumerate(top):
        entry = dict(cfg)
        entry['status'] = 'running'
        t0 = time.time()
        # build hnsw with given parameters
        try:
            knn = hnsw.build_knn(profiles, k=k, ef_construction=cfg['ef_construction'], M=cfg['M'], ef_search=cfg['ef_search'])
            t1 = time.time()
            entry['build_time_actual'] = t1 - t0
            if knn.nnz == 0:
                entry['error'] = 'empty knn'
                entry['status'] = 'failed'
            else:
                # compute recall between exact and approx
                approx_neighbors = [list(knn.tocsr().indices[knn.tocsr().indptr[i]:knn.tocsr().indptr[i+1]]) for i in range(n)]
                # recall as fraction of exact neighbors retrieved
                total = 0; hit = 0
                for i in range(n):
                    ref = set(exact_knn[i])
                    appr = set(approx_neighbors[i][:k])
                    total += len(ref)
                    hit += len(ref & appr)
                entry['recall'] = hit / total if total>0 else 0.0

                # rerank and compute MST edge overlap
                reranked = rerank_candidates_with_exact(profiles, knn)
                approx_edges = mst_edge_set_from_csr(reranked)
                entry['mst_edge_overlap'] = len(exact_edges & approx_edges) / len(exact_edges) if len(exact_edges)>0 else 0.0
                entry['status'] = 'done'
        except Exception as exc:
            entry['status'] = 'error'
            entry['error'] = str(exc)
        results['configs'].append(entry)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(json.dumps(results, indent=2))
    print('Wrote', out)

if __name__ == '__main__':
    main()

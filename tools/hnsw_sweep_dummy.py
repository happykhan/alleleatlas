"""HNSW parameter sweep on synthetic no-missing cgMLST-like data.

This script runs multiple synthetic datasets and sweeps HNSW parameters
(ef_search, ef_construction, M) and hashed encoding dims, collecting
recall@k, precision@k, MRR@k, MST overlap, and timings. It also tries
unioning HNSW and MinHash candidates and measures rerank costs.

Outputs JSON to out/hnsw_sweep.json
"""
import json
import time
from pathlib import Path
import numpy as np
import hashlib
import math

from alleleatlas import distances_exact
from alleleatlas.backends import hnsw, minhash_lsh
import scipy.sparse as sp
import tools.testbed as tb

OUT = Path('out')
OUT.mkdir(exist_ok=True)

# experiment parameters
N_DATASETS = 10
n = 100
d = 50
K = 5
seeds = list(range(100, 100 + N_DATASETS))

M_vals = [16, 32]
ef_construction_vals = [100, 200]
ef_search_vals = [50, 100, 200, 500]
hash_dims = [4096, 8192, 16384]
candidate_multipliers = [1, 3, 5, 10]

results = []


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
                # fallback to blake2b
                h = hashlib.blake2b(key.encode('utf8')).hexdigest()
            idx = int(h, 16) % dim
            X[i, idx] = 1.0
    return X


for seed in seeds:
    np.random.seed(seed)
    profiles = np.random.randint(1, 10, size=(n, d)).astype(float)

    # compute exact baseline
    legacy = distances_exact.get_legacy_pairwise(profiles, method='numba')
    D = legacy[:, :, 0].astype(float)
    exact_knn = distances_exact.knn_from_profiles(profiles, k=K)
    exact_edges = tb.mst_edge_set_from_csr(exact_knn)

    for dim in hash_dims:
        Xh = one_hot_hash_encode(profiles, dim=dim, hash_name='md5')
        for Mval in M_vals:
            for efc in ef_construction_vals:
                for efs in ef_search_vals:
                    for mult in candidate_multipliers:
                        candidate_k = min(n - 1, mult * K)

                        # build HNSW on hashed encoding with candidate_k
                        t0 = time.time()
                        try:
                            knn_h = hnsw.build_knn(Xh, candidate_k, ef_construction=efc, M=Mval, ef_search=efs)
                        except Exception as exc:
                            results.append({'seed': seed, 'dim': dim, 'M': Mval, 'efc': efc, 'efs': efs, 'mult': mult, 'error': str(exc)})
                            continue
                        t1 = time.time()
                        build_time = t1 - t0

                        # extract candidate lists
                        if hasattr(knn_h, 'tocsr'):
                            approx_neighbors = [list(knn_h.tocsr().indices[knn_h.tocsr().indptr[i]:knn_h.tocsr().indptr[i+1]]) for i in range(knn_h.shape[0])]
                        else:
                            approx_neighbors = knn_h

                        # rerank candidates to final K using exact distances and time rerank
                        t_r0 = time.time()
                        nrows = profiles.shape[0]
                        rows, cols, data = [], [], []
                        for i in range(nrows):
                            cand = approx_neighbors[i][:candidate_k]
                            if len(cand) == 0:
                                continue
                            vec = profiles[i]
                            dists = np.array([distances_exact._pair_distance_row(vec, profiles)[j] for j in cand], dtype=float)
                            # pick top-K after rerank
                            if len(dists) > K:
                                sel = np.argpartition(dists, K)[:K]
                                chosen = [cand[s] for s in sel]
                                chosen_d = dists[sel]
                            else:
                                chosen = cand
                                chosen_d = dists
                            for j, dd in zip(chosen, chosen_d):
                                rows.append(i)
                                cols.append(int(j))
                                data.append(float(dd))
                        t_r1 = time.time()
                        rerank_time = t_r1 - t_r0
                        if len(rows) == 0:
                            reranked_csr = sp.csr_matrix((n, n), dtype=float)
                        else:
                            reranked_csr = sp.csr_matrix((data, (rows, cols)), shape=(n, n))

                        # compute metrics
                        recall = tb.recall_between_knn(exact_knn, approx_neighbors, k=K)
                        precision = None
                        try:
                            precision = sum([len(set(approx_neighbors[i][:K]) & set(exact_knn.tocsr().indices[exact_knn.tocsr().indptr[i]:exact_knn.tocsr().indptr[i+1]]))/K for i in range(n)])/n
                        except Exception:
                            precision = None
                        # MRR
                        def mrr(ref, approx, k=K):
                            refcsr = ref.tocsr()
                            n = refcsr.shape[0]
                            total = 0.0
                            for i in range(n):
                                refset = set(refcsr.indices[refcsr.indptr[i]:refcsr.indptr[i+1]])
                                rank_score = 0.0
                                for r, c in enumerate(approx[i][:k], start=1):
                                    if c in refset:
                                        rank_score = 1.0 / r
                                        break
                                total += rank_score
                            return total / n

                        mrrv = mrr(exact_knn, approx_neighbors, k=K)

                        # MST overlap after rerank
                        try:
                            approx_edges = tb.mst_edge_set_from_csr(reranked_csr)
                            overlap = len(exact_edges & approx_edges) / len(exact_edges) if len(exact_edges) > 0 else 0.0
                        except Exception:
                            overlap = None

                        results.append({
                            'seed': seed,
                            'dim': dim,
                            'M': Mval,
                            'ef_construction': efc,
                            'ef_search': efs,
                            'candidate_multiplier': mult,
                            'candidate_k': candidate_k,
                            'build_time': build_time,
                            'rerank_time': rerank_time,
                            'recall_at_k': recall,
                            'precision_at_k': precision,
                            'mrr_at_k': mrrv,
                            'mst_edge_overlap': overlap,
                            'nnz_candidates': int(knn_h.nnz) if hasattr(knn_h, 'nnz') else None,
                        })

                        # now union with MinHash candidates and rerank
                        try:
                            t0 = time.time()
                            knn_m = minhash_lsh.build_knn(profiles, candidate_k, num_perm=256)
                            t1 = time.time()
                            mh_time = t1 - t0
                            if hasattr(knn_m, 'tocsr'):
                                approx_neighbors_m = [list(knn_m.tocsr().indices[knn_m.tocsr().indptr[i]:knn_m.tocsr().indptr[i+1]]) for i in range(knn_m.shape[0])]
                            else:
                                approx_neighbors_m = knn_m

                            # union per-query
                            t_r0 = time.time()
                            rows_u, cols_u, data_u = [], [], []
                            approx_union = []
                            for i in range(nrows):
                                union_c = list(dict.fromkeys(approx_neighbors[i] + approx_neighbors_m[i]))[:candidate_k]
                                approx_union.append(union_c)
                                if len(union_c) == 0:
                                    continue
                                vec = profiles[i]
                                dists = np.array([distances_exact._pair_distance_row(vec, profiles)[j] for j in union_c], dtype=float)
                                if len(dists) > K:
                                    sel = np.argpartition(dists, K)[:K]
                                    chosen = [union_c[s] for s in sel]
                                    chosen_d = dists[sel]
                                else:
                                    chosen = union_c
                                    chosen_d = dists
                                for j, dd in zip(chosen, chosen_d):
                                    rows_u.append(i)
                                    cols_u.append(int(j))
                                    data_u.append(float(dd))
                            t_r1 = time.time()
                            rerank_time_u = t_r1 - t_r0
                            if len(rows_u) == 0:
                                reranked_csr_u = sp.csr_matrix((n, n), dtype=float)
                            else:
                                reranked_csr_u = sp.csr_matrix((data_u, (rows_u, cols_u)), shape=(n, n))

                            recall_u = tb.recall_between_knn(exact_knn, approx_union, k=K)
                            mrr_u = mrr(exact_knn, approx_union, k=K)
                            try:
                                approx_edges_u = tb.mst_edge_set_from_csr(reranked_csr_u)
                                overlap_u = len(exact_edges & approx_edges_u) / len(exact_edges) if len(exact_edges) > 0 else 0.0
                            except Exception:
                                overlap_u = None

                            results.append({
                                'seed': seed,
                                'dim': dim,
                                'M': Mval,
                                'ef_construction': efc,
                                'ef_search': efs,
                                'candidate_multiplier': mult,
                                'candidate_k': candidate_k,
                                'method': 'hnsw+minhash_union',
                                'build_time': build_time,
                                'minhash_time': mh_time,
                                'rerank_time': rerank_time_u,
                                'recall_at_k': recall_u,
                                'mrr_at_k': mrr_u,
                                'mst_edge_overlap': overlap_u,
                            })

                        except Exception as exc:
                            results.append({'seed': seed, 'dim': dim, 'M': Mval, 'efc': efc, 'efs': efs, 'mult': mult, 'union_error': str(exc)})

# write results
out_file = OUT / 'hnsw_sweep.json'
out_file.write_text(json.dumps(results, indent=2), encoding='utf-8')
print('Wrote', out_file)

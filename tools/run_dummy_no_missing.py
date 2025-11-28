"""Generate a 100x50 no-missing synthetic profile and run the testbed logic.

This script mirrors `tools.testbed` but uses an in-memory generated dataset
and prints a compact JSON summary to stdout and writes `out/dummy_no_missing.json`.
"""
import json
import time
from pathlib import Path
import numpy as np

from alleleatlas import distances_exact
from alleleatlas.backends import hamming, minhash_lsh, sklearn_nn, hnsw, rerank
import scipy.sparse as sp

# import helper functions from testbed
import tools.testbed as tb
import hashlib

OUT_DIR = Path('out')
OUT_DIR.mkdir(exist_ok=True)

np.random.seed(42)
# generate alleles in 1..9 (no zeros to avoid legacy missing semantics)
n = 100
d = 50
profiles = np.random.randint(1, 10, size=(n, d)).astype(float)
K = 5

results = {'input_rows': n, 'input_cols': d, 'runs': {}}

# exact (legacy)
start = time.time()
legacy_arr = distances_exact.get_legacy_pairwise(profiles, method='numba')
end = time.time()
D = legacy_arr[:, :, 0].astype(float)
results['runs']['exact'] = {'time_s': end - start, 'shape': D.shape}

# exact kNN
start = time.time()
exact_knn = distances_exact.knn_from_profiles(profiles, k=K)
end = time.time()
results['runs']['exact_knn'] = {'time_s': end - start}

backends = [
    ('hamming', hamming),
    ('minhash_lsh', minhash_lsh),
    ('sklearn_nn', sklearn_nn),
    ('hnsw', hnsw),
    ('hnsw_rerank', rerank),
]

for name, mod in backends:
    try:
        if hasattr(mod, 'build_knn'):
            t0 = time.time()
            knn = mod.build_knn(profiles, K)
            t1 = time.time()
            results['runs'][name] = {'time_s': t1 - t0}
            # compute recall vs exact
            if hasattr(knn, 'tocsr'):
                approx_neighbors = [list(knn.tocsr().indices[knn.tocsr().indptr[i]:knn.tocsr().indptr[i+1]]) for i in range(knn.shape[0])]
            else:
                approx_neighbors = knn
            recall = tb.recall_between_knn(exact_knn, approx_neighbors, k=K)
            results['runs'][name]['recall'] = recall

            # rerank candidates with exact distance
            if hasattr(knn, 'tocsr'):
                approx_csr = knn.tocsr()
                n = approx_csr.shape[0]
                rows, cols, data = [], [], []
                for i in range(n):
                    cand = approx_csr.indices[approx_csr.indptr[i]: approx_csr.indptr[i+1]]
                    if cand.size == 0:
                        continue
                    vec = profiles[i]
                    dists = np.array([distances_exact._pair_distance_row(vec, profiles)[j] for j in cand], dtype=float)
                    for j, d in zip(cand, dists):
                        rows.append(i)
                        cols.append(int(j))
                        data.append(float(d))
                if len(rows) == 0:
                    reranked_csr = sp.csr_matrix((n, n), dtype=float)
                else:
                    reranked_csr = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
                try:
                    exact_edges = tb.mst_edge_set_from_csr(exact_knn)
                    approx_edges = tb.mst_edge_set_from_csr(reranked_csr)
                    overlap = len(exact_edges & approx_edges) / len(exact_edges) if len(exact_edges) > 0 else 0.0
                except Exception:
                    overlap = None
                results['runs'][name]['mst_edge_overlap'] = overlap
        elif name == 'hnsw_rerank' and hasattr(mod, 'hnsw_rerank'):
            t0 = time.time()
            knn = mod.hnsw_rerank(profiles, K)
            t1 = time.time()
            results['runs'][name] = {'time_s': t1 - t0}
            approx_neighbors = [list(knn.tocsr().indices[knn.tocsr().indptr[i]:knn.tocsr().indptr[i+1]]) for i in range(knn.shape[0])]
            recall = tb.recall_between_knn(exact_knn, approx_neighbors, k=K)
            results['runs'][name]['recall'] = recall
            approx_edges = tb.mst_edge_set_from_csr(knn)
            exact_edges = tb.mst_edge_set_from_csr(exact_knn)
            overlap = len(exact_edges & approx_edges) / len(exact_edges) if len(exact_edges) > 0 else 0.0
            results['runs'][name]['mst_edge_overlap'] = overlap
    except Exception as exc:
        results['runs'][name] = {'skipped': True, 'reason': str(exc)}


def precision_at_k(reference_csr, approx_neighbors, k=5):
    ref = reference_csr.tocsr()
    n = ref.shape[0]
    precisions = []
    for i in range(n):
        ref_set = set(ref.indices[ref.indptr[i]: ref.indptr[i+1]])
        appr = approx_neighbors[i][:k]
        if len(appr) == 0:
            precisions.append(0.0)
            continue
        hit = len(set(appr) & ref_set)
        precisions.append(hit / len(appr))
    return float(np.mean(precisions))


def mrr_at_k(reference_csr, approx_neighbors, k=5):
    ref = reference_csr.tocsr()
    n = ref.shape[0]
    rr = []
    for i in range(n):
        ref_set = set(ref.indices[ref.indptr[i]: ref.indptr[i+1]])
        found = 0.0
        for rank, cand in enumerate(approx_neighbors[i][:k], start=1):
            if cand in ref_set:
                found = 1.0 / rank
                break
        rr.append(found)
    return float(np.mean(rr))


def one_hot_hash_encode(profiles, dim=2048):
    """Encode integer profiles (n, d) into binary hashed one-hot vectors (n, dim).

    Each (locus, value) is hashed to a bucket in [0, dim) using md5 for
    deterministic hashing, and that bucket is set to 1. Missing values are
    ignored.
    """
    n, d = profiles.shape
    X = np.zeros((n, dim), dtype=float)
    for i in range(n):
        for j in range(d):
            v = profiles[i, j]
            # skip NaN
            if isinstance(v, float) and np.isnan(v):
                continue
            key = f"{j}_{int(v)}"
            h = hashlib.md5(key.encode('utf8')).hexdigest()
            idx = int(h, 16) % dim
            X[i, idx] = 1.0
    return X


# Additional HNSW experiment: raw vs hashed encoding
if 'hnsw' in results['runs'] and not results['runs']['hnsw'].get('skipped'):
    try:
        # raw already captured as 'hnsw'
        # now run HNSW on one-hot hashed encoding
        dim = 4096
        Xh = one_hot_hash_encode(profiles, dim=dim)
        t0 = time.time()
        knn_h = hnsw.build_knn(Xh, K)
        t1 = time.time()
        key = 'hnsw_onehot'
        results['runs'][key] = {'time_s': t1 - t0}
        if hasattr(knn_h, 'tocsr'):
            approx_neighbors_h = [list(knn_h.tocsr().indices[knn_h.tocsr().indptr[i]:knn_h.tocsr().indptr[i+1]]) for i in range(knn_h.shape[0])]
        else:
            approx_neighbors_h = knn_h
        recall_h = tb.recall_between_knn(exact_knn, approx_neighbors_h, k=K)
        prec_h = precision_at_k(exact_knn, approx_neighbors_h, k=K)
        mrr_h = mrr_at_k(exact_knn, approx_neighbors_h, k=K)
        results['runs'][key]['recall'] = recall_h
        results['runs'][key]['precision_at_k'] = prec_h
        results['runs'][key]['mrr_at_k'] = mrr_h
        # rerank and MST overlap
        if hasattr(knn_h, 'tocsr'):
            approx_csr = knn_h.tocsr()
            n = approx_csr.shape[0]
            rows, cols, data = [], [], []
            for i in range(n):
                cand = approx_csr.indices[approx_csr.indptr[i]: approx_csr.indptr[i+1]]
                if cand.size == 0:
                    continue
                vec = profiles[i]
                dists = np.array([distances_exact._pair_distance_row(vec, profiles)[j] for j in cand], dtype=float)
                for j, d in zip(cand, dists):
                    rows.append(i)
                    cols.append(int(j))
                    data.append(float(d))
            if len(rows) == 0:
                reranked_csr = sp.csr_matrix((n, n), dtype=float)
            else:
                reranked_csr = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
            try:
                exact_edges = tb.mst_edge_set_from_csr(exact_knn)
                approx_edges = tb.mst_edge_set_from_csr(reranked_csr)
                overlap = len(exact_edges & approx_edges) / len(exact_edges) if len(exact_edges) > 0 else 0.0
            except Exception:
                overlap = None
            results['runs'][key]['mst_edge_overlap'] = overlap
    except Exception as exc:
        results['runs']['hnsw_onehot'] = {'skipped': True, 'reason': str(exc)}

# write JSON
out_path = OUT_DIR / 'dummy_no_missing.json'
out_path.write_text(json.dumps(results, indent=2), encoding='utf-8')
print('Wrote', out_path)
print(json.dumps(results, indent=2))

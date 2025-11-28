#!/usr/bin/env python3
"""Compute reranked recall@k for existing out/*.json runs and update files.

This script looks for `out/hybrid_union_verify.json`, `out/faiss_binary_verify.json`,
and `out/hnsw_onehot_verify.json`. For each run entry it recomputes final reranked
top-k lists by requesting candidate lists (or rebuilding with stored cfg) and
computes `recall_reranked` and writes the updated JSON back.
"""
from pathlib import Path
import json
import numpy as np
import scipy.sparse as sp
from alleleatlas import distances_exact


def load_profiles(profile='test_trees/paraC.profile'):
    from alleleatlas.convert_to_parquet import load_profiles
    df, _ = load_profiles(profile)
    if 'ST' in df.columns:
        df = df.drop(columns=['ST'])
    return df.values.astype(int)


def rerank_and_compute(profiles, approx_csr, k=5):
    n = approx_csr.shape[0]
    rows, cols, data = [], [], []
    for i in range(n):
        cand = approx_csr.indices[approx_csr.indptr[i]: approx_csr.indptr[i+1]]
        if cand.size == 0:
            continue
        vec = profiles[i]
        dists = np.array([distances_exact._pair_distance_row(vec, profiles)[j] for j in cand], dtype=float)
        if len(dists) > k:
            sel = np.argpartition(dists, k)[:k]
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
        return 0.0, sp.csr_matrix((profiles.shape[0], profiles.shape[0]), dtype=float)
    csr = sp.csr_matrix((data, (rows, cols)), shape=(profiles.shape[0], profiles.shape[0]), dtype=float)
    # compute recall: count how many of exact top-k are present in reranked neighbors
    legacy = distances_exact.get_legacy_pairwise(profiles, method='numba')
    exact_D = legacy[:, :, 0].astype(float)
    exact_knn = np.argsort(exact_D, axis=1)[:, 1:6]
    n = profiles.shape[0]
    total = 0
    hit = 0
    rows2 = csr.tocsr()
    for i in range(n):
        ref = set(exact_knn[i])
        appr = set(rows2.indices[rows2.indptr[i]: rows2.indptr[i+1]])
        total += len(ref)
        hit += len(ref & appr)
    recall_reranked = hit/total if total>0 else 0.0
    return recall_reranked, csr


def update_file(path: Path):
    j = json.loads(path.read_text())
    if 'input_rows' not in j:
        return
    profiles = load_profiles()
    if path.name == 'hnsw_onehot_verify.json':
        # These runs were done without storing reranked lists; we need to rebuild per-cfg
        from alleleatlas.backends import hnsw
        runs = j['runs']
        for r in runs:
            cfg = r.get('cfg', {})
            dim = cfg.get('dim', 4096)
            # build hashed encoding
            from tools.run_onehot_hnsw_on_profile import one_hot_hash_encode
            Xh = one_hot_hash_encode(profiles, dim=dim)
            k_cand = min(profiles.shape[0]-1, cfg.get('candidate_k', 5))
            knn_h = hnsw.build_knn(Xh, k_cand, ef_construction=cfg.get('ef_construction',100), M=cfg.get('M',16), ef_search=cfg.get('ef_search',100))
            recall_reranked, csr = rerank_and_compute(profiles, knn_h, k=5)
            r['recall_reranked'] = recall_reranked
    else:
        # For hybrid and faiss we expect candidate CSR already produced by their runners
        runs = j.get('runs', [])
        for r in runs:
            # we need to reconstruct the candidate CSR; the runners wrote files but didn't persist candidate lists
            # try to re-run the same runner logic via available backends
            cfg = r.get('cfg', {})
            # attempt to use hybrid backend
            from alleleatlas.backends import hybrid_union, faiss_binary
            if path.name == 'hybrid_union_verify.json':
                approx = hybrid_union.build_knn(profiles, k=5, minhash_perms=cfg.get('minhash_perms',256), minhash_k=cfg.get('minhash_k',50), faiss_dim_bits=cfg.get('faiss_dim_bits',16384), use_faiss_ivf=cfg.get('use_faiss_ivf',False))
            elif path.name == 'faiss_binary_verify.json':
                approx = faiss_binary.build_knn(profiles, k=5, dim_bits=cfg.get('dim_bits',16384), use_ivf=cfg.get('use_ivf',False))
            else:
                approx = None
            if approx is None:
                r['recall_reranked'] = None
            else:
                recall_reranked, csr = rerank_and_compute(profiles, approx, k=5)
                r['recall_reranked'] = recall_reranked
    path.write_text(json.dumps(j, indent=2))
    print('Updated', path)


def main():
    out = Path('out')
    for name in ['hybrid_union_verify.json', 'faiss_binary_verify.json', 'hnsw_onehot_verify.json']:
        p = out / name
        if p.exists():
            update_file(p)


if __name__ == '__main__':
    main()

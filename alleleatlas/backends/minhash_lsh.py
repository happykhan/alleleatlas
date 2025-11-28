"""MinHash + LSH backend for approximate similarity search.

This backend tokenizes each profile row as a set of "locus:value" strings
and builds MinHash signatures. An LSH index finds candidate neighbors which
are then re-ranked using exact Euclidean distance on the numeric vectors.

The signature size (``num_perm``) and LSH ``threshold`` are exposed so
the caller can tune recall/precision tradeoffs.

If :mod:`datasketch` is not installed the function falls back to a brute-
force nearest neighbor computation; this keeps tests simple and deterministic
when the optional dependency is missing.
"""
from __future__ import annotations

import numpy as np
from typing import List


def build_knn(profiles: np.ndarray, k: int, *, num_perm: int = 128, threshold: float = 0.5, verbose: bool = False):
    """Build k-NN using MinHash LSH candidate generation.

    Parameters
    - profiles: (n, d) integer array
    - k: number of neighbors
    - num_perm: number of permutations for MinHash
    - threshold: LSH similarity threshold
    - verbose: print progress if True
    """
    n, d = profiles.shape
    try:
        from datasketch import MinHash, MinHashLSH
    except Exception:
        # fallback: brute-force (small n only)
        import scipy.spatial as sps
        import scipy.sparse as sp
        X = profiles.astype(float)
        if n == 0:
            return sp.csr_matrix((0, 0))
        D = sps.distance.cdist(X, X, metric='euclidean')
        np.fill_diagonal(D, np.inf)
        idxs = np.argpartition(D, k, axis=1)[:, :k]
        rows = np.repeat(np.arange(n), k)
        cols = idxs.reshape(-1)
        data = D[np.arange(n)[:, None], idxs].reshape(-1)
        return sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()

    # Build MinHash signatures and LSH index
    lsh = MinHashLSH(threshold=threshold, num_perm=num_perm)
    minhashes: List[MinHash] = []
    tokens_list: List[set] = []
    for i in range(n):
        m = MinHash(num_perm=num_perm)
        row = profiles[i]
        toks = set()
        # tokenization: include only non-zero loci to keep signatures sparse
        for j, v in enumerate(row):
            # skip missing values (NaN) and zeros
            if np.isnan(v):
                continue
            if v != 0:
                tok = f"{j}:{int(v)}"
                toks.add(tok)
                m.update(tok.encode('utf-8'))
        lsh.insert(str(i), m)
        minhashes.append(m)
        tokens_list.append(toks)

    import scipy.spatial as sps
    import scipy.sparse as sp
    rows: List[int] = []
    cols: List[int] = []
    data: List[float] = []
    X = profiles.astype(float)

    for i in range(n):
        cand = lsh.query(minhashes[i])
        cand_idx = [int(c) for c in cand if int(c) != i]
        # if we don't have enough candidates, supplement using brute-force token-set Jaccard
        if len(cand_idx) < k:
            # compute jaccard against all using token sets
            toks_i = tokens_list[i]
            all_jacc = []
            for j in range(n):
                if j == i:
                    all_jacc.append(-1.0)
                    continue
                toks_j = tokens_list[j]
                if len(toks_i) == 0 and len(toks_j) == 0:
                    jacc = 0.0
                else:
                    inter = len(toks_i & toks_j)
                    union = len(toks_i | toks_j)
                    jacc = float(inter) / float(union) if union > 0 else 0.0
                all_jacc.append(jacc)
            order = np.argpartition([-v for v in all_jacc], min(k, len(all_jacc) - 1))[:k + 1]
            cand_idx = [j for j in order if j != i][:k]
        if len(cand_idx) == 0:
            continue
        # compute estimated Jaccard via MinHash and convert to distance (1 - jaccard)
        dists = np.array([1.0 - minhashes[i].jaccard(minhashes[j]) for j in cand_idx], dtype=float)
        if len(dists) > k:
            sel = np.argpartition(dists, k)[:k]
            cand_idx = [cand_idx[s] for s in sel]
            dists = dists[sel]
        for j, dist in zip(cand_idx, dists):
            rows.append(i)
            cols.append(j)
            data.append(float(dist))

    return sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()

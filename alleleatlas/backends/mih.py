"""Simple MIH-like exact retrieval helper for small Hamming radii.

This is a minimal fallback implementation intended for moderate-size datasets
and small Hamming radii. It is not a production-grade MIH but provides a
fast path for low-radius exact neighbor retrieval by brute-force on packed bits
for moderate N.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def pack_bitset_from_profiles(profiles: np.ndarray, dim_bits: int = 16384):
    # reuse simple packing similar to faiss backend
    n, d = profiles.shape
    nbytes = dim_bits // 8
    out = np.zeros((n, nbytes), dtype='uint8')
    import hashlib
    hfn = hashlib.blake2b
    for i in range(n):
        for j, v in enumerate(profiles[i]):
            try:
                if isinstance(v, float) and np.isnan(v):
                    continue
            except Exception:
                pass
            key = f"{j}_{int(v)}".encode('utf8')
            h = hfn(key).digest()
            idx = int.from_bytes(h, 'big') % dim_bits
            out[i, idx // 8] |= (1 << (idx % 8))
    return out


def hamming_popcount(u: np.ndarray, v: np.ndarray) -> int:
    # u and v are uint8 arrays; compute popcount of xor
    xor = np.bitwise_xor(u, v)
    # vectorized popcount via table
    table = np.unpackbits(np.arange(256, dtype='uint8')[:, None], axis=1).sum(axis=1)
    return table[xor].sum()


def build_knn_exact_by_radius(profiles: np.ndarray, k: int, dim_bits: int = 16384, radius: int = 10):
    """Return CSR of neighbors within Hamming radius (exact brute-force on packed bits).

    For moderate N (<= 100k) and small radius this can be fast.
    """
    packed = pack_bitset_from_profiles(profiles, dim_bits=dim_bits)
    n = packed.shape[0]
    rows = []
    cols = []
    data = []
    # simple O(n^2) loop - only suitable for small-to-moderate n
    for i in range(n):
        ui = packed[i]
        for j in range(n):
            if i == j:
                continue
            dist = hamming_popcount(ui, packed[j])
            if dist <= radius:
                rows.append(i)
                cols.append(j)
                data.append(float(dist))
    if len(rows) == 0:
        return sp.csr_matrix((n, n), dtype=float)
    return sp.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=float)

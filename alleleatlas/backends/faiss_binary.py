"""Faiss binary backend for allele profiles.

This backend encodes allele (col,value) pairs into a fixed-length bit vector via hashing
and builds a FAISS binary index (BinaryFlat or BinaryIVF). It returns a CSR matrix of
approximate neighbors in the same shape as other backends.

This module is optional: if `faiss` is not installed the backend raises ImportError with
an actionable message.
"""
from __future__ import annotations

import hashlib
import math
from typing import Optional

import numpy as np
import scipy.sparse as sp

try:
    import faiss
except Exception:
    faiss = None


def _onebit_hash_pack(profiles: np.ndarray, dim_bits: int = 16384, hash_name: str = 'blake2b') -> np.ndarray:
    """Encode integer profile array (n, d) to packed binary array (n, dim_bits/8) of dtype uint8.

    Uses a stable hash of "col_value" to set a single bit for that allele state.
    Missing values (np.nan) are skipped.
    """
    n, d = profiles.shape
    if dim_bits % 8 != 0:
        raise ValueError('dim_bits must be multiple of 8')
    nbytes = dim_bits // 8
    out = np.zeros((n, nbytes), dtype='uint8')
    hashfn = hashlib.blake2b if hash_name == 'blake2b' else hashlib.md5
    for i in range(n):
        row = profiles[i]
        for j, v in enumerate(row):
            if v is None:
                continue
            # allow NaN floats
            try:
                if isinstance(v, float) and np.isnan(v):
                    continue
            except Exception:
                pass
            key = f"{j}_{int(v)}".encode('utf8')
            h = hashfn(key).digest()
            # interpret hash as big int
            idx = int.from_bytes(h, 'big') % dim_bits
            byte = idx // 8
            bit = idx % 8
            out[i, byte] |= (1 << bit)
    return out


def _packed_to_faiss_bytes(packed_uint8: np.ndarray) -> np.ndarray:
    # faiss expects contiguous uint8 array of shape (n, nb)
    return np.ascontiguousarray(packed_uint8)


def build_knn(profiles: np.ndarray, k: int, *, dim_bits: int = 16384, use_ivf: bool = False,
              nlist: Optional[int] = None, nprobe: int = 1, hash_name: str = 'blake2b') -> sp.csr_matrix:
    """Build k-NN candidate lists using Faiss binary index.

    Returns a CSR matrix where each row i contains candidate neighbor indices (not reranked).
    If `faiss` is not installed raises ImportError with instructions.
    """
    if faiss is None:
        raise ImportError('faiss is not available. Install with `pip install faiss-cpu` or use conda faiss packages.')

    n, d = profiles.shape
    if n == 0:
        return sp.csr_matrix((0, 0), dtype=float)

    packed = _onebit_hash_pack(profiles, dim_bits=dim_bits, hash_name=hash_name)
    xb = _packed_to_faiss_bytes(packed)
    # ensure dtype and contiguity
    xb = np.ascontiguousarray(xb, dtype='uint8')
    expected_nbytes = dim_bits // 8
    if xb.shape[1] != expected_nbytes:
        raise ValueError(f'packed byte width {xb.shape[1]} does not match expected {expected_nbytes} for dim_bits={dim_bits}')

    # build index
    if not use_ivf:
        index = faiss.IndexBinaryFlat(dim_bits)
        index.add(xb)
        # sanity check
        if index.ntotal != n:
            # something went wrong adding vectors
            raise RuntimeError(f'faiss index ntotal={index.ntotal} expected {n}')
        # search
        D, Ind = index.search(xb, min(k, n - 1))
    else:
        # compute default nlist
        if nlist is None:
            nlist = max(32, int(math.sqrt(n)))
        quantizer = faiss.IndexBinaryFlat(dim_bits)
        index = faiss.IndexBinaryIVF(quantizer, dim_bits, nlist)
        index.nprobe = nprobe
        index.train(xb)
        index.add(xb)
        if index.ntotal != n:
            raise RuntimeError(f'faiss index ntotal={index.ntotal} expected {n}')
        D, Ind = index.search(xb, min(k, n - 1))

    # FAISS uses -1 for missing neighbors; convert to CSR
    rows = []
    cols = []
    data = []
    for i in range(n):
        # some FAISS builds may return invalid indices in edge cases; clamp them
        row_inds = [int(x) for x in Ind[i]]
        for j in row_inds:
            if j < 0 or j >= n:
                continue
            if j == i:
                continue
            rows.append(i)
            cols.append(int(j))
            data.append(1.0)
    if len(rows) == 0:
        return sp.csr_matrix((n, n), dtype=float)
    return sp.csr_matrix((data, (rows, cols)), shape=(n, n), dtype=float)

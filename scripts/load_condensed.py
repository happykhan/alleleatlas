#!/usr/bin/env python3
"""Load condensed distance memmaps created by compute_distances_condensed.py

Provides helpers to:
 - read the meta JSON produced by the batch jobs
 - expand condensed memmap to a full symmetric matrix (for small N)
 - build a k-NN sparse adjacency (CSR) from condensed memmap without full expansion

Usage examples:
  # Expand to dense (only if N small / you have RAM)
  python3 scripts/load_condensed.py --meta tmp/st313.meta.json --mode full --out tmp/st313.full.npy

  # Build k-NN sparse neighbor graph and save as .npz
  python3 scripts/load_condensed.py --meta tmp/st313.meta.json --mode knn --k 50 --out tmp/st313.knn.npz

"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import math
from scipy.sparse import coo_matrix, csr_matrix, save_npz


def read_meta(meta_path: Path):
    with open(meta_path, 'r') as fh:
        meta = json.load(fh)
    return meta


def condensed_index(n: int, i: int, j: int) -> int:
    # index into condensed array for i<j
    return n * i - (i * (i + 1) // 2) + (j - i - 1)


def load_condensed_memmap(path: str):
    # load memmap as 1D int32 array
    arr = np.memmap(path, dtype=np.int32, mode='r')
    return arr


def condensed_to_full(dist_mem: np.memmap, n: int) -> np.ndarray:
    """Expand condensed 1D upper-triangle array to full symmetric nxn matrix.

    Warning: requires O(n^2) RAM.
    """
    if n * n > 5000 * 5000:
        raise MemoryError(f"Refusing to expand matrix with n={n} (would allocate {n*n} entries). Use knn mode instead.")
    mat = np.zeros((n, n), dtype=np.float32)
    # fill upper triangle
    for i in range(n):
        for j in range(i + 1, n):
            idx = condensed_index(n, i, j)
            mat[i, j] = float(dist_mem[idx])
            mat[j, i] = mat[i, j]
    return mat


def build_knn_sparse(dist_mem: np.memmap, n: int, k: int) -> csr_matrix:
    """Build a k-NN sparse matrix (symmetric) from condensed distances without expanding fully.

    This computes, for each row i, its k nearest neighbors by scanning distances.
    Complexity: O(n^2) reads but memory efficient.
    """
    rows = []
    cols = []
    data = []
    for i in range(n):
        # gather distances to all j != i
        # build index array for j>i and j<i
        dist_row = np.empty((n,), dtype=np.float32)
        dist_row[:] = np.inf
        # j < i : condensed index (j,i)
        if i > 0:
            js = np.arange(0, i)
            idxs = (n * js - (js * (js + 1) // 2) + (i - js - 1)).astype(np.int64)
            dist_row[js] = dist_mem[idxs]
        # j > i : condensed index (i,j)
        if i < n - 1:
            js = np.arange(i + 1, n)
            idxs = (n * i - (i * (i + 1) // 2) + (js - i - 1)).astype(np.int64)
            dist_row[js] = dist_mem[idxs]

        # find k nearest (exclude inf)
        k_use = min(k, n - 1)
        nbr_idx = np.argpartition(dist_row, k_use)[:k_use]
        for j in nbr_idx:
            if not np.isfinite(dist_row[j]):
                continue
            rows.append(i)
            cols.append(j)
            data.append(float(dist_row[j]))

    # make symmetric by adding transpose
    rows_sym = np.array(rows + cols, dtype=np.int32)
    cols_sym = np.array(cols + rows, dtype=np.int32)
    data_sym = np.array(data + data, dtype=np.float32)
    coo = coo_matrix((data_sym, (rows_sym, cols_sym)), shape=(n, n))
    return coo.tocsr()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--meta', type=Path, required=True, help='Path to meta JSON produced by batch job')
    parser.add_argument('--mode', choices=['full', 'knn'], default='knn')
    parser.add_argument('--k', type=int, default=50, help='k for k-NN (knn mode)')
    parser.add_argument('--mutual', action='store_true', help='Keep only mutual (reciprocal) k-NN edges')
    parser.add_argument('--out', type=Path, required=True, help='Output file (.npy for full, .npz for knn)')
    args = parser.parse_args()

    meta = read_meta(args.meta)
    dist0 = meta.get('dist0')
    if dist0 is None:
        raise RuntimeError('meta does not contain dist0')
    names = meta.get('names')

    # infer n from names if present, else from memmap length
    if names and Path(names).exists():
        with open(names, 'r') as fh:
            n = sum(1 for _ in fh)
    else:
        # infer from len of memmap: mlen = n*(n-1)//2
        mm = np.memmap(dist0, dtype=np.int32, mode='r')
        mlen = mm.shape[0]
        # solve n*(n-1)//2 = mlen
        n = int((1 + math.sqrt(1 + 8 * mlen)) / 2)

    dist_mem = load_condensed_memmap(dist0)

    if args.mode == 'full':
        print(f'Expanding condensed memmap to full matrix (n={n})...')
        mat = condensed_to_full(dist_mem, n)
        np.save(str(args.out), mat)
        print(f'Wrote dense matrix to {args.out}')
    else:
        print(f'Building k-NN (k={args.k}) sparse matrix (n={n})...')
        csr = build_knn_sparse(dist_mem, n, args.k)
        if args.mutual:
            # keep only mutual edges: (i,j) where both csr[i,j] and csr[j,i] exist
            csr = csr.tocsr()
            csr_T = csr.transpose().tocsr()
            # elementwise minimum keeps only reciprocal edges (non-zero in both)
            mutual = csr.minimum(csr_T)
            csr = mutual.tocsr()
        save_npz(str(args.out), csr)
        print(f'Wrote k-NN sparse matrix to {args.out}')


if __name__ == '__main__':
    main()

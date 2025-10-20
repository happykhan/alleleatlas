#!/usr/bin/env python3
"""Compute condensed distance arrays (upper-triangle) from a cgMLST profile.

This script implements a memory-efficient, chunked distance computation:
 - writes mat, mask and counts to disk as memmaps
 - creates condensed memmaps for distances (two channels)
 - uses a Pool(initializer=...) so workers attach memmaps once
 - numba-parallel kernels compute per-row blocks which are written to condensed

Usage:
  python3 scripts/compute_distances_condensed.py INPUT_PROFILE out_prefix --nproc 4

Output:
  out_prefix.dist0.memmap  (int32, condensed upper-triangle)
  out_prefix.dist1.memmap  (int32, condensed upper-triangle)
  out_prefix.mat.memmap    (int16) -- optional
  out_prefix.names.txt     (sample names)

This is intended as a drop-in alternative to the SharedArray-based getDistance
implementation and is safe to run on systems with limited RAM (uses disk-backed memmaps).
"""

from __future__ import annotations
import argparse
import math
import os
import sys
import tempfile
import numpy as np
import multiprocessing as mp
from multiprocessing import Pool
from pathlib import Path
from typing import Tuple
from numba import njit, prange

# When running this script directly (or under SLURM) the package
# `alleleatlas` may not be on PYTHONPATH. Add the repository root to
# sys.path so local imports work even without PYTHONPATH set.
try:
    from alleleatlas.cluster.pHierCC import prepare_mat
except ModuleNotFoundError:
    repo_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(repo_root))
    from alleleatlas.cluster.pHierCC import prepare_mat


def condensed_index(n: int, i: int, j: int) -> int:
    """Index into condensed array for i < j.

    idx = n*i - i*(i+1)//2 + (j - i - 1)
    """
    return n * i - (i * (i + 1) // 2) + (j - i - 1)


@njit(parallel=True, fastmath=True)
def dual_dist_kernel(mat, present_mask, present_count, s, e, allowed_missing):
    N = mat.shape[0]
    L = mat.shape[1]
    out0 = np.zeros((e - s, N), dtype=np.int32)
    out1 = np.zeros((e - s, N), dtype=np.int32)
    am = np.float64(allowed_missing)
    Lf = np.float64(L)
    eps = np.float64(1e-12)
    for ii in prange(s, e):
        i_idx = ii - s
        ql = np.float64(present_count[ii])
        for j in range(ii):
            rl = np.float64(present_count[j])
            ad = np.float64(0.0)
            al = eps
            for k in range(L):
                if present_mask[j, k]:
                    if present_mask[ii, k]:
                        al += np.float64(1.0)
                        if mat[ii, k] != mat[j, k]:
                            ad += np.float64(1.0)
            ll2 = ql - am * Lf
            if ll2 > al:
                ad += (ll2 - al)
                al = ll2
            out1[i_idx, j] = np.int32(ad / al * Lf + np.float64(0.5))
            ll = max(ql, rl) - am * Lf
            if ll > al:
                ad += (ll - al)
                al = ll
            out0[i_idx, j] = np.int32(ad / al * Lf + np.float64(0.5))
    return out0, out1


@njit(parallel=True, fastmath=True)
def p_dist_kernel(mat, present_mask, present_count, s, e):
    N = mat.shape[0]
    L = mat.shape[1]
    out = np.zeros((e - s, N), dtype=np.int32)
    Lf = np.float64(L)
    for ii in prange(s, e):
        i_idx = ii - s
        for j in range(ii):
            ad = np.float64(0.0)
            al = np.float64(0.0)
            for k in range(L):
                if present_mask[j, k] and present_mask[ii, k]:
                    al += np.float64(1.0)
                    if mat[ii, k] != mat[j, k]:
                        ad += np.float64(1.0)
            denom = al + np.float64(1.0)
            tmp = (ad + np.float64(0.5)) / denom
            if tmp >= np.float64(1.0):
                tmp = np.float64(1.0) - np.float64(1e-12)
            val = -np.log(np.float64(1.0) - tmp) * Lf * np.float64(100.0) + np.float64(0.5)
            out[i_idx, j] = np.int32(val)
    return out


# Worker-global state (memmaps attached in initializer)
_STATE = {
    'mat': None,
    'mask': None,
    'count': None,
    'dist0': None,
    'dist1': None,
    'n': None,
}


def worker_init(mat_path, mat_shape, mat_dtype, mask_path, mask_shape, mask_dtype,
                count_path, count_shape, count_dtype, dist0_path, dist1_path, n):
    """Attach memmaps once per worker."""
    global _STATE
    _STATE['mat'] = np.memmap(mat_path, dtype=np.dtype(mat_dtype), mode='r', shape=tuple(mat_shape))
    _STATE['mask'] = np.memmap(mask_path, dtype=np.dtype(mask_dtype), mode='r', shape=tuple(mask_shape))
    _STATE['count'] = np.memmap(count_path, dtype=np.dtype(count_dtype), mode='r', shape=tuple(count_shape))
    mlen = (n * (n - 1)) // 2
    _STATE['dist0'] = np.memmap(dist0_path, dtype=np.int32, mode='r+', shape=(mlen,))
    _STATE['dist1'] = np.memmap(dist1_path, dtype=np.int32, mode='r+', shape=(mlen,))
    _STATE['n'] = n


def worker_task(args: Tuple[str, int, int, float]):
    """Compute rows [s,e) and write j<i distances into condensed memmaps.

    args = (func, s, e, allowed_missing)
    """
    func, s, e, allowed_missing = args
    mat = _STATE['mat']
    mask = _STATE['mask']
    count = _STATE['count']
    n = _STATE['n']
    dist0 = _STATE['dist0']
    dist1 = _STATE['dist1']

    if func == 'dual':
        out0, out1 = dual_dist_kernel(mat, mask, count, s, e, allowed_missing)
    else:
        out0 = p_dist_kernel(mat, mask, count, s, e)
        out1 = np.zeros_like(out0)

    # write condensed: for each ii in [s,e), for j in 0..ii-1 write index(j,ii)
    for ii in range(s, e):
        i_idx = ii - s
        for j in range(ii):
            idx = condensed_index(n, j, ii)
            dist0[idx] = int(out0[i_idx, j])
            dist1[idx] = int(out1[i_idx, j])
    return (s, e)


def compute_condensed(profile_path: str, out_prefix: str, nproc: int = 4, allowed_missing: float = 0.05):
    console = None
    profile_path = str(profile_path)
    # prepare matrix and names
    mat, names = prepare_mat(profile_path)
    n = mat.shape[0]
    L = mat.shape[1] - 1 if mat.shape[1] > 1 else 0

    # mat from prepare_mat includes id column; keep full mat for algorithm compatibility
    # convert to compact dtype
    if mat.dtype != np.int16:
        mat = mat.astype(np.int16, copy=False)

    present_mask = (mat > 0).astype(np.uint8)
    present_count = np.sum(mat > 0, axis=1).astype(np.int32)

    # write memmaps
    out_prefix = str(out_prefix)
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    mat_path = out_prefix + '.mat.memmap'
    mask_path = out_prefix + '.mask.memmap'
    count_path = out_prefix + '.count.memmap'
    names_path = out_prefix + '.names.txt'

    # create memmaps and write
    mp_mat = np.memmap(mat_path, dtype=mat.dtype, mode='w+', shape=mat.shape)
    mp_mat[:] = mat[:]
    mp_mask = np.memmap(mask_path, dtype=present_mask.dtype, mode='w+', shape=present_mask.shape)
    mp_mask[:] = present_mask[:]
    mp_count = np.memmap(count_path, dtype=present_count.dtype, mode='w+', shape=present_count.shape)
    mp_count[:] = present_count[:]
    with open(names_path, 'w') as fh:
        for s in names:
            fh.write(str(s) + '\n')

    # condensed lengths
    mlen = (n * (n - 1)) // 2
    dist0_path = out_prefix + '.dist0.memmap'
    dist1_path = out_prefix + '.dist1.memmap'
    mp.dist0 = np.memmap(dist0_path, dtype=np.int32, mode='w+', shape=(mlen,))
    mp.dist0[:] = 0
    mp.dist1 = np.memmap(dist1_path, dtype=np.int32, mode='w+', shape=(mlen,))
    mp.dist1[:] = 0

    # create chunks (one per worker)
    n_workers = max(1, min(nproc, mp.cpu_count()))
    base = n // n_workers
    extras = n % n_workers
    ranges = []
    s = 0
    for i in range(n_workers):
        e = s + base + (1 if i < extras else 0)
        if e > s:
            ranges.append((s, e))
        s = e

    # prepare init args
    init_args = (mat_path, mat.shape, mat.dtype.str,
                 mask_path, present_mask.shape, present_mask.dtype.str,
                 count_path, present_count.shape, present_count.dtype.str,
                 dist0_path, dist1_path, n)

    # run pool
    pool = Pool(processes=n_workers, initializer=worker_init, initargs=init_args)
    try:
        tasks = [('dual', s, e, allowed_missing) for s, e in ranges]
        for _ in pool.imap_unordered(worker_task, tasks):
            pass
    finally:
        pool.close()
        pool.join()

    print(f"Wrote condensed memmaps: {dist0_path}, {dist1_path} (n={n})")
    return dist0_path, dist1_path, mat_path, names_path


def main():
    parser = argparse.ArgumentParser(description="Compute condensed distances (memmap) from cgMLST profile")
    parser.add_argument('profile')
    parser.add_argument('out_prefix')
    parser.add_argument('--nproc', type=int, default=4)
    parser.add_argument('--allowed-missing', type=float, default=0.05)
    args = parser.parse_args()

    compute_condensed(args.profile, args.out_prefix, nproc=args.nproc, allowed_missing=args.allowed_missing)


if __name__ == '__main__':
    main()

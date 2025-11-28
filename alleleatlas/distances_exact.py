"""Exact cgMLST distance utilities (refactored from legacy getDistance).

This module provides behaviorally-identical operations for exact
pairwise distance computation between cgMLST profiles while exposing a
cleaner API. Missing data (NaN) is ignored when comparing loci.

Functions:
- pairwise_distances(profiles: np.ndarray) -> np.ndarray
- knn_from_profiles(profiles: np.ndarray, k: int) -> scipy.sparse.csr_matrix

The distance metric counts the number of differing allele calls across
positions where both profiles have non-missing values.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from typing import Optional
import numba as nb


@nb.jit(nopython=True)
def legacy_dual_numba(mat, s, e, allowed_missing=0.05):
    """Port of legacy dual_dist from original getDistance implementation.

    Returns an int32 array shaped (e-s, mat.shape[0], 2) where [:,:,0] and
    [:,:,1] follow the legacy semantics.
    """
    dist = np.zeros((e - s, mat.shape[0], 2), dtype=np.int32)
    n_loci = mat.shape[1]
    for i in range(s, e):
        ql = 0
        for kk in range(mat.shape[1]):
            if mat[i, kk] > 0:
                ql += 1
        for j in range(i):
            rl = 0.0
            ad = 1e-4
            al = 1e-4
            for k in range(n_loci):
                if mat[j, k] > 0:
                    rl += 1.0
                    if mat[i, k] > 0:
                        al += 1.0
                        if mat[i, k] != mat[j, k]:
                            ad += 1.0
            ll = max(ql, rl) - allowed_missing * n_loci
            ll2 = ql - allowed_missing * n_loci

            if ll2 > al:
                ad += ll2 - al
                al = ll2
            # dist[i-s, j, 1] = int(ad/al * n_loci + 0.5)
            dist[i - s, j, 1] = np.int32(ad / al * n_loci + 0.5)

            if ll > al:
                ad += ll - al
                al = ll
            dist[i - s, j, 0] = np.int32(ad / al * n_loci + 0.5)
    return dist


# legacy p_dist not needed for current API; removed to keep file focused


def pairwise_legacy_numba(profiles: np.ndarray, allowed_missing: float = 0.05) -> np.ndarray:
    """Compute the legacy-style mismatch/overlap matrix by calling the
    ported numba `legacy_dual_dist` and expanding the triangular output
    into a full symmetric (n, n, 2) array.

    This function treats missing values as zeros (<=0) to match the
    original implementation semantics.
    """
    # prepare integer matrix where missing (NaN) -> 0
    P = profiles.astype(float)
    # treat NaN as 0 (missing)
    P[np.isnan(P)] = 0.0
    M = P.astype(np.int32)
    n = M.shape[0]
    # call the numba port; it returns shape (n, n, 2) with filled lower-triangle
    tri = legacy_dual_numba(M, 0, n, allowed_missing)
    # tri has shape (n, n, 2) where entries for j<i are set in tri[i-s, j]
    out = np.zeros((n, n, 2), dtype=np.int32)
    # fill lower triangle and mirror to upper triangle
    for i in range(n):
        for j in range(i):
            out[i, j, 0] = tri[i, j, 0]
            out[i, j, 1] = tri[i, j, 1]
            out[j, i, 0] = tri[i, j, 0]
            out[j, i, 1] = tri[i, j, 1]
    # diagonal: mismatch 0, comparable loci = number of positive entries in row
    for i in range(n):
        out[i, i, 0] = 0
        out[i, i, 1] = int((M[i] > 0).sum())
    return out


def pairwise_legacy_python(profiles: np.ndarray, allowed_missing: float = 0.05) -> np.ndarray:
    """Pure-Python implementation of the legacy dual_dist semantics.

    Returns (n, n, 2) int32 array where entries match the numba port
    produced by `legacy_dual_dist` expanded by `pairwise_using_legacy`.
    """
    P = profiles.astype(float)
    P[np.isnan(P)] = 0.0
    M = P.astype(np.int32)
    n, n_loci = M.shape
    out = np.zeros((n, n, 2), dtype=np.int32)
    for i in range(n):
        ql = int((M[i] > 0).sum())
        for j in range(i):
            rl = 0
            ad = 1e-4
            al = 1e-4
            for k in range(n_loci):
                if M[j, k] > 0:
                    rl += 1
                    if M[i, k] > 0:
                        al += 1
                        if M[i, k] != M[j, k]:
                            ad += 1
            ll = max(ql, rl) - allowed_missing * n_loci
            ll2 = ql - allowed_missing * n_loci

            if ll2 > al:
                ad += ll2 - al
                al = ll2
            v1 = int(ad / al * n_loci + 0.5)

            if ll > al:
                ad += ll - al
                al = ll
            v0 = int(ad / al * n_loci + 0.5)

            out[i, j, 0] = v0
            out[i, j, 1] = v1
            out[j, i, 0] = v0
            out[j, i, 1] = v1
    for i in range(n):
        out[i, i, 0] = 0
        out[i, i, 1] = int((M[i] > 0).sum())
    return out


def get_legacy_pairwise(profiles: np.ndarray, method: str = 'numba', allowed_missing: float = 0.05) -> np.ndarray:
    """Dispatch to the legacy implementation.

    method: 'numba' to use the numba port, 'python' to use the pure-python
    reimplementation. Raises ValueError for unknown methods.
    """
    if method == 'numba':
        return pairwise_legacy_numba(profiles, allowed_missing=allowed_missing)
    elif method == 'python':
        return pairwise_legacy_python(profiles, allowed_missing=allowed_missing)
    else:
        raise ValueError(f'unknown legacy method: {method}')


def _pair_distance_row(x: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Compute distances from x to every row in Y.

    x: (d,), Y: (n, d) -> returns (n,) distances (int counts). Missing
    values (np.nan) are ignored per-position.
    """
    # mask of valid positions for each comparison: (~np.isnan(x)) & (~np.isnan(Y))
    valid = ~np.isnan(Y) & (~np.isnan(x))
    # number of comparable loci per row
    counts = valid.sum(axis=1)
    # compute mismatches where both valid and values differ
    # Use broadcasting: compare Y to x where valid
    diffs = (Y != x) & valid
    mismatches = diffs.sum(axis=1)
    # For pairs with zero comparable loci, set distance to np.inf
    mismatches = mismatches.astype(float)
    mismatches[counts == 0] = np.nan
    return mismatches


# NOTE: legacy behavior is authoritative for "exact" distances.
# The previous NaN-aware `pairwise_distances` implementation was removed
# to avoid divergent semantics. Use `get_legacy_pairwise(..., method=...)`
# for exact pairwise matrices. The helper `_pair_distance_row` remains for
# internal use but is NOT the canonical "exact" implementation.


def knn_from_profiles(profiles: np.ndarray, k: int, method: str = 'numba', allowed_missing: float = 0.05) -> sp.csr_matrix:
    """Build a k-NN sparse distance matrix using the legacy pairwise metric.

    By default this calls the numba legacy port via
    `get_legacy_pairwise(profiles, method='numba')`. Use `method='python'`
    to select the pure-Python legacy reimplementation.

    Returns a square (n, n) CSR matrix where entry (i, j) is the distance
    from i to its neighbor j for the k nearest neighbors of i. Self-neighbors
    are not included.
    """
    # obtain legacy pairwise (n, n, 2) where [:,:,0] is the distance
    legacy = get_legacy_pairwise(profiles, method=method, allowed_missing=allowed_missing)
    D = legacy[:, :, 0].astype(float)
    n = D.shape[0]
    if n == 0:
        return sp.csr_matrix((0, 0), dtype=float)
    # mask self to avoid selecting
    np.fill_diagonal(D, np.inf)
    k_eff = max(0, min(k, n - 1))
    if k_eff == 0:
        return sp.csr_matrix((n, n), dtype=float)
    # use argpartition per-row for performance
    idxs = np.argpartition(D, k_eff, axis=1)[:, :k_eff]
    # gather values
    rr = np.repeat(np.arange(n), k_eff)
    cc = idxs.reshape(-1)
    vals = D[np.arange(n)[:, None], idxs].reshape(-1)
    mat = sp.csr_matrix((vals, (rr, cc)), shape=(n, n))
    return mat


def main_cli(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description='Exact cgMLST distances')
    parser.add_argument('parquet', help='input parquet with normalized profiles')
    parser.add_argument('--k', type=int, default=50)
    parser.add_argument('--out', help='output .npz for k-NN', required=True)
    args = parser.parse_args(argv)
    import pandas as pd
    df = pd.read_parquet(args.parquet)
    if 'ST' in df.columns:
        loci = df.drop(columns=['ST'])
    else:
        loci = df.iloc[:, 1:]
    X = loci.to_numpy(dtype=float)
    import time as _time
    t0 = _time.time()
    knn = knn_from_profiles(X, args.k)
    t1 = _time.time()
    import scipy.sparse as sp
    sp.save_npz(args.out, knn)
    print(f'Wrote {args.out} (n={X.shape[0]},d={X.shape[1]}) k={args.k} time={t1-t0:.2f}s')


if __name__ == '__main__':
    main_cli()


def pairwise_mismatch_and_overlap(profiles: np.ndarray, *, missing_value: Optional[int] = None) -> np.ndarray:
    """Return an (n, n, 2) integer array where [:,:,0] is mismatch counts
    and [:,:,1] is number of comparable loci between pairs.

    If `missing_value` is provided (e.g., 0), entries equal to that value
    are treated as missing. If `missing_value` is None, NaN is treated as
    missing (the expected internal representation used elsewhere).
    """
    n, d = profiles.shape
    out = np.zeros((n, n, 2), dtype=np.int32)
    # prepare float view with NaNs for missing if needed
    if missing_value is not None:
        # treat missing_value (e.g., 0) as missing -> convert to NaN
        P = profiles.astype(float)
        P[P == missing_value] = np.nan
    else:
        P = profiles.astype(float)

    for i in range(n):
        valid = ~np.isnan(P) & (~np.isnan(P[i]))
        counts = valid.sum(axis=1).astype(np.int32)
        diffs = ((P != P[i]) & valid).sum(axis=1).astype(np.int32)
        out[i, :, 0] = diffs
        out[i, :, 1] = counts
    # zero diagonal: no mismatches and d comparable (or 0 if missing)
    for i in range(n):
        out[i, i, 0] = 0
        # out[i,i,1] remains counts (will typically be number of non-missing loci)
    return out

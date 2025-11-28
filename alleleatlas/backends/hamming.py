"""Hamming-distance backend for allele profiles.

Builds k-NN where distance = count of differing loci (ignoring NaNs).
Returns a scipy.sparse.csr_matrix with distances for each row to its k nearest neighbors.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def build_knn(profiles: np.ndarray, k: int):
    """Compute Hamming-like distances and return k-NN sparse matrix.

    profiles: (n, d) float array with NaN for missing entries.
    """
    n, d = profiles.shape
    if n == 0:
        return sp.csr_matrix((0, 0))
    rows = []
    cols = []
    data = []
    # compute pairwise distance row-by-row to avoid O(n^2) memory in dense
    for i in range(n):
        xi = profiles[i]
        # compute differing counts vs all rows: treat NaN as missing and ignore those loci
        valid_mask = ~np.isnan(xi)
        # for rows with all NaNs, set distance to large / skip
        if not np.any(valid_mask):
            # no valid loci to compare; skip
            continue
        # compute differences only on valid_mask columns
        sub = profiles[:, valid_mask]
        # create boolean array: differs if not equal (and not NaN in other)
        # treat NaN in j as missing -> ignore those positions
        other_notnan = ~np.isnan(sub)
        # compute per-row valid positions: both not nan
        both_valid = other_notnan[:, :]  # (n, m)
        xi_sub = xi[valid_mask]
        # broadcast compare
        diffs = (sub != xi_sub) & both_valid
        # count differing loci per row and compute overlap (number of compared loci)
        counts = np.asarray(diffs.sum(axis=1)).astype(float)
        overlap = np.asarray(both_valid.sum(axis=1)).astype(float)
        # compute normalized distance = fraction differing (0..1); if no overlap -> inf
        with np.errstate(divide='ignore', invalid='ignore'):
            frac = np.where(overlap > 0, counts / overlap, np.inf)
        # ignore self
        frac[i] = np.inf
        # select up to k smallest by fractional distance
        m = min(k, max(1, len(frac) - 1))
        if m <= 0:
            continue
        idxs = np.argpartition(frac, m)[:m]
        sel = idxs[np.argsort(frac[idxs])]
        for j in sel:
            if np.isinf(frac[j]):
                continue
            rows.append(i)
            cols.append(int(j))
            data.append(float(frac[j]))
    if len(rows) == 0:
        return sp.csr_matrix((n, n))
    M = sp.coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
    return M

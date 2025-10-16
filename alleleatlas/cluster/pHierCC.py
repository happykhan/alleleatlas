#!/usr/bin/env python

# pHierCC.py
# pipeline for Hierarchical Clustering of cgMLST
#
# Author: Zhemin Zhou
# Lisence: GPLv3
#
# New assignment: pHierCC -p <allelic_profile> -o <output_prefix>
# Incremental assignment: pHierCC -p <allelic_profile> -o <output_prefix> -i <old_cluster_npz>
# Input format (tab delimited):
# ST_id gene1 gene2
# 1 1 1
# 2 1 2
# ...

import gzip
import pandas as pd
import numpy as np
from pathlib import Path
from multiprocessing import Pool #, set_start_method
from scipy.spatial import distance as ssd
from scipy.cluster.hierarchy import linkage
from rich.console import Console

from alleleatlas.cluster.getDistance import getDistance

console = Console()

def prepare_mat(profile_file, use_parquet_cache=True) :
    """Read an allelic profile file and prepare a numeric matrix for clustering.

    Parameters
    ----------
    profile_file : str or path-like
        Path to a tab-delimited allelic profile file. Expected format is a
        header row (optional) followed by rows like: ST_id allele1 allele2 ...
        The first column is treated as the sequence/sample identifier.
    
    use_parquet_cache : bool, optional
        If True, check for .parquet cache file and use it if available for faster
        loading. If not available, load from TSV and create cache. Default True.

    Returns
    -------
    mat : ndarray of int, shape (n_samples, n_loci+1)
        Numeric matrix where the first column is an integer id (either the
        original ST id if numeric, or an assigned sequential index) and the
        remaining columns are integer-coded alleles. Missing/invalid alleles
        are converted to 0.

    names : ndarray
        The original first-column values (as strings) that can be used to
        map rows back to sample identifiers.

    Notes
    -----
    - If allele values are already integers this function uses them directly
      and drops rows with non-positive ST ids.
    - If allele values are non-numeric strings the function maps unique
      allele tokens to integers (deterministically via np.unique). Empty
      strings, '0', or tokens starting with '-' are treated as missing and
      mapped to 0.
    """
    
    # Try parquet cache first for faster loading
    if use_parquet_cache:
        parquet_file = str(profile_file) + '.parquet'
        try:
            if Path(parquet_file).exists():
                console.print(f'[green]✓[/green] Loading cached parquet: {parquet_file}')
                df = pd.read_parquet(parquet_file)
                names = df.iloc[:, 0].values
                mat = df.values.astype(int)
                return mat, names
        except Exception as e:
            console.print(f'[yellow]Warning:[/yellow] Could not load parquet cache: {e}')

    # Read file as strings to tolerate mixed allele formats and preserve header
    mat = pd.read_csv(profile_file, sep='\t', header=None, dtype=str, na_filter=False).values

    # Determine which columns to keep: always keep the first column (ids),
    # and drop any columns in the header row that start with '#'.
    allele_columns = np.array([i == 0 or (not h.startswith('#')) for i, h in enumerate(mat[0])])

    # Remove header row and any unwanted comment columns
    mat = mat[1:, allele_columns]

    try :
        # Fast path: all fields are numeric -> cast and filter out non-positive ids
        mat = mat.astype(int)
        mat = mat[mat.T[0] > 0]
        names = mat.T[0]
    except Exception:
        # Fallback for non-numeric allele tokens:
        # - preserve original names
        # - replace first column with sequential integer ids (1..n)
        names = mat.T[0].copy()
        mat.T[0] = np.arange(1, mat.shape[0]+1)

        # Map allele string tokens -> integers.
        # - Treat empty, '0', or tokens starting with '-' as missing (0).
        # - Preserve purely numeric tokens by converting them to their int values.
        # - Assign new positive integer codes to non-numeric tokens, starting
        #   after the maximum numeric allele value to avoid collisions.
        alleles = mat[:, 1:]
        uniques = np.unique(alleles)

        # Collect numeric tokens (as ints) and non-numeric tokens separately.
        # Guard against extremely large numeric-looking strings (they may be
        # hashes or concatenated fields) by only accepting numeric tokens
        # within a reasonable range.
        numeric_vals = {}
        non_numeric = []
        MAX_NUM = 10 ** 9
        for tag in uniques:
            if tag in {'', '0'} or tag.startswith('-'):
                continue
            try:
                iv = int(tag)
                # accept numeric tokens only if within expected bounds
                if 0 <= iv <= MAX_NUM:
                    numeric_vals[tag] = iv
                else:
                    non_numeric.append(tag)
            except (ValueError, OverflowError):
                non_numeric.append(tag)

        # Determine the starting code for non-numeric tokens. Start after the
        # largest numeric allele value (if any), otherwise start at 1. 0 is
        # reserved for missing values.
        next_code = (max(numeric_vals.values()) + 1) if numeric_vals else 1

        # Build mapping dict deterministically using the order from uniques.
        d = {}
        for tag in uniques:
            if tag in {'', '0'} or tag.startswith('-'):
                d[tag] = 0
            elif tag in numeric_vals:
                d[tag] = numeric_vals[tag]
            else:
                d[tag] = next_code
                next_code += 1

        # Apply mapping safely with an explicit loop to avoid numpy overflow
        # when converting large Python ints to C longs.
        rows, cols = mat.shape[0], mat.shape[1] - 1
        mapped = np.zeros((rows, cols), dtype=int)
        for i in range(rows):
            for j in range(cols):
                mapped[i, j] = d.get(mat[i, j+1], 0)

        # Reconstruct mat: first column is the sequential ids we set earlier
        mat = np.hstack((mat[:, [0]].astype(int), mapped))

    # Clamp any negative values to zero (missing)
    mat[mat < 0] = 0
    
    # Save to parquet cache for faster future loading
    if use_parquet_cache:
        try:
            parquet_file = str(profile_file) + '.parquet'
            df_cache = pd.DataFrame(mat)
            df_cache.to_parquet(parquet_file, compression='snappy', index=False)
            console.print(f'[green]✓[/green] Cached profile to parquet: {parquet_file}')
        except Exception as e:
            console.print(f'[yellow]Warning:[/yellow] Could not save parquet cache: {e}')
    
    return mat, names

def phierCC(profile, output, append, n_proc, allowed_missing, precomputed_dist=None):
    '''pHierCC takes a file containing allelic profiles (as in https://pubmlst.org/data/) and works
    out hierarchical clusters of the full dataset based on a minimum-spanning tree.
    
    Parameters:
        profile: path to profile file
        output: output prefix
        append: append mode (file path or empty string)
        n_proc: number of processes
        allowed_missing: allowed missing value threshold
        precomputed_dist: optional pre-computed distance matrix (3D array from getDistance)
    '''
    pool = Pool(n_proc)

    profile_file, cluster_file, old_cluster = profile, output + '.npz', append

    mat, names = prepare_mat(profile_file)
    n_loci = mat.shape[1] - 1

    console.log(f'Loaded in allelic profiles with dimension: {mat.shape[0]} and {mat.shape[1]}. The first column is assumed to be type id.')
    console.log('Start HierCC assignments')

    # prepare existing clusters
    # Initialize cls to satisfy static analyzers; it will be set if append is True
    cls = None

    if not append:
        absence = np.sum(mat <= 0, 1)
        mat[:] = mat[np.argsort(absence, kind='mergesort')]
        typed = {}
    else :
        od = np.load(old_cluster, allow_pickle=True)
        cls = od['hierCC']
        try :
            n = od['names']
        except Exception:
            n = cls.T[0]
        typed = {c: id for id, c in enumerate(n)}
    if len(typed) > 0:
        console.log(f'Loaded in {len(typed)} old HierCC assignments.')
        # mat_idx = np.array([t in typed for t in names])
        mat_idx = np.argsort([typed.get(t, len(typed)) for t in names])
        mat[:] = mat[mat_idx]
        names[:] = names[mat_idx]
        start = np.sum([t in typed for t in names])
        if names.dtype != np.int64 :
            mat.T[0] = np.arange(1, mat.shape[0]+1)
    else :
        start = 0

    res = np.repeat(mat.T[0], int(mat.shape[1]) + 1).reshape(mat.shape[0], -1)
    res[res < 0] = np.max(mat.T[0]) + 100
    res.T[0] = mat.T[0]
    
    # Use precomputed distance if provided, otherwise compute
    if precomputed_dist is not None:
        console.log('Using precomputed distance matrix')
        dist = precomputed_dist
    else:
        console.log('Calculate distance matrix')
        dist = getDistance(mat, 'dual_dist', pool, start, allowed_missing)
    if append :
        for n, r in zip(names, res) :
            if n in typed and cls is not None:
                r[:] = cls[typed[n]]
    else :
        dist[:, :, 0] += dist[:, :, 0].T
        console.log('Start Single linkage clustering')
        slc = linkage(ssd.squareform(dist[:, :, 0]), method='single')
        index = { s:i for i, s in enumerate(mat.T[0]) }
        # Use empty lists for future internal nodes (instead of None) to avoid
        # static typing issues while still allowing concatenation when nodes
        # are created.
        descendents = [ [m] for m in mat.T[0] ] + [[] for _ in np.arange(mat.shape[0]-1)]
        for idx, c in enumerate(slc.astype(int)) :
            n_id = idx + mat.shape[0]
            # choose the child with the lexicographically smaller first element
            d = sorted([int(c[0]), int(c[1])], key=lambda x: descendents[x][0])
            min_id = min(descendents[d[0]])
            descendents[n_id] = descendents[d[0]] + descendents[d[1]]
            for tgt in descendents[d[1]] :
                res[index[tgt], c[2]+1:] = res[index[min_id], c[2]+1:]

    console.log('Attach genomes onto the tree.')
    for id, (r, d) in enumerate(zip(res[start:], dist[:, :, 1])):
        if id + start > 0 :
            i = np.argmin(d[:id+start])
            min_d = d[i]
            if r[min_d + 1] > res[i, min_d + 1]:
                r[min_d + 1:] = res[i, min_d + 1:]
    res.T[0] = mat.T[0]
    res = res[np.argsort(res.T[0])]
    np.savez_compressed(cluster_file, hierCC=res, names=names)

    with gzip.open(output + '.HierCC.gz', 'wt') as fout:
        fout.write('#ST_id\t{0}\n'.format('\t'.join(['HC' + str(id) for id in np.arange(n_loci+1)])))
        for n, r in zip(names, res):
            fout.write('\t'.join([str(n)] + [str(rr) for rr in r[1:]]) + '\n')

    console.log(f'NPZ  clustering result (for production mode): {output}.npz')
    console.log(f'TEXT clustering result (for visual inspection and HCCeval): {output}.HierCC.gz')
    pool.close()


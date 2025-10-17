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

import numpy as np
from rich.console import Console
from sklearn.metrics import silhouette_score, normalized_mutual_info_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from multiprocessing import Pool
import SharedArray as sa
from tempfile import NamedTemporaryFile
from alleleatlas.cluster.getDistance import getDistance

console = Console()


def prepare_mat(profile_file) :
    """Read an allelic profile file and prepare a numeric matrix for clustering.

    Parameters
    ----------
    profile_file : str or path-like
        Path to a tab-delimited allelic profile file. Expected format is a
        header row (optional) followed by rows like: ST_id allele1 allele2 ...
        The first column is treated as the sequence/sample identifier.

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
    return mat, names


def get_similarity2(data) :
    method, cc1, cc2 = data
    if np.unique(cc1).size == 1 and  np.unique(cc1).size == 1 :
        return 1.
    return method(cc1, cc2)

def get_similarity(method, cluster, stepwise, pool) :
    console.log('Calculating NMIs...')
    similarity = np.ones([cluster.shape[1], cluster.shape[1]], dtype=np.float64)
    for i1, cc1 in enumerate(cluster.T) :
        if i1 % 10 == 0 :
            console.log(f'    NMIs between level {i1 * stepwise} and greater levels')
        similarity[i1, i1+1:] = pool.map(get_similarity2, [ [method, cc1, cc2] for cc2 in cluster.T[i1+1:] ])
        similarity[i1+1:, i1] = similarity[i1, i1+1:]
    similarity[similarity>0.999] = 0.999
    similarity[similarity<0.0] = 0.0
    return similarity

def get_silhouette(profile, cluster, stepwise, pool, precomputed_dist=None) :
    if precomputed_dist is not None:
        console.log('Using precomputed distance matrix')
        dist = precomputed_dist
    else:
        console.log('Calculating pairwise distance ...')
        dist = getDistance(profile, 'p_dist', pool)
    with NamedTemporaryFile(dir='.', prefix='HCCeval_') as file :
        dist_buf = 'file://{0}.dist'.format(file.name)
        dist2 = sa.create(dist_buf, dist.shape[:2], dist.dtype)
        dist2[:] = dist[:, :, 0] + dist[:, :, 0].T
        del dist
        console.log('Calculating Silhouette score ...')
        silhouette = np.array(pool.map(get_silhouette2, [ [dist_buf, tag] for tag in cluster.T ]))
        sa.delete(dist_buf)
    return silhouette

def get_silhouette2(data) :
    dist_buf, tag = data
    s = np.unique(tag)
    if 2 <= s.size < tag.shape[0] :
        dist = sa.attach(dist_buf)
        ss = silhouette_score(dist.astype(float), tag, metric = 'precomputed')
        return ss
    else :
        return 0.


def my_phierCC(mat, names, distance, output, n_proc):
    '''pHierCC takes a file containing allelic profiles (as in https://pubmlst.org/data/) and works
    out hierarchical clusters of the full dataset based on a minimum-spanning tree.
    
    Parameters:

    '''
    pool = Pool(n_proc)
    # Prepare the allele profile matrix
    n_loci = mat.shape[1] - 1
    cluster_file = str(output) + '.npz'
    console.log(f'Loaded in allelic profiles with dimension: {mat.shape[0]} and {mat.shape[1]}. The first column is assumed to be type id.')
    console.log('Start HierCC assignments')

    absence = np.sum(mat <= 0, 1)
    mat[:] = mat[np.argsort(absence, kind='mergesort')]
    start = 0

    res = np.repeat(mat.T[0], int(mat.shape[1]) + 1).reshape(mat.shape[0], -1)
    res[res < 0] = np.max(mat.T[0]) + 100
    res.T[0] = mat.T[0]
    dist = distance
    
    # ============================================================================
    # SINGLE LINKAGE HIERARCHICAL CLUSTERING
    # ============================================================================
    # This section builds a hierarchical dendrogram tree from the distance matrix
    # using single linkage (minimum distance criterion). The result assigns each
    # sample to a sequence of hierarchical clusters at increasing distance thresholds.
    
    # STEP 1: Symmetrize the distance matrix
    # The distance matrix dist[:, :, 0] is lower-triangular. We need a full 
    # symmetric matrix for linkage analysis. Add its transpose to make it symmetric.
    dist[:, :, 0] += dist[:, :, 0].T
    console.log('Start Single linkage clustering')
    
    # STEP 2: Compute hierarchical clustering linkage tree
    # ssd.squareform() converts lower-triangle to condensed distance vector.
    # scipy's linkage() builds a dendrogram using single linkage (minimum distance
    # between clusters). Output: (n_samples-1) x 4 array where each row is:
    #   [cluster1_idx, cluster2_idx, merge_distance, sample_count]
    slc = linkage(ssd.squareform(dist[:, :, 0]), method='single')
    
    # STEP 3: Create reverse mapping: sample_id -> row_index
    # Example: if mat.T[0] = [1001, 1002, 1003], then index = {1001:0, 1002:1, 1003:2}
    # This allows fast lookup of sample position in the result matrix.
    index = { s:i for i, s in enumerate(mat.T[0]) }
    
    # STEP 4: Initialize descendant tracking lists
    # For each internal node in the dendrogram tree, track which original samples
    # (leaves) are its descendants. Start with each sample as a leaf in its own list.
    # Then add empty lists for future internal nodes (will be filled as tree builds).
    # descendents[i] for i < n_samples = list of leaf samples under that leaf
    # descendents[i] for i >= n_samples = list of leaf samples under that internal node
    descendents = [ [m] for m in mat.T[0] ] + [[] for _ in np.arange(mat.shape[0]-1)]
    
    # STEP 5: Traverse linkage tree and propagate cluster assignments
    # For each merge event in the dendrogram (from bottom to top):
    
    for idx, c in enumerate(slc.astype(int)) :
        n_id = idx + mat.shape[0]
        d = sorted([int(c[0]), int(c[1])], key=lambda x:descendents[x][0])
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
    return res



def evalHCC(mat, names, cluster, output, stepwise, n_proc, precomputed_dist=None) :
    '''evalHCC evaluates a HierCC scheme using varied statistic summaries.
    
    Parameters:
        profile: profile matrix
        cluster: cluster matrix
        output: output prefix
        stepwise: stepwise parameter
        n_proc: number of processes
        precomputed_dist: optional pre-computed distance matrix (3D array from getDistance)
    '''
    '''evalHCC evaluates a HierCC scheme using varied statistic summaries.'''
    pool = Pool(n_proc)

    # Profile should be names + profile matrix
    profile = np.column_stack((names, mat))
    cluster = np.column_stack((names, cluster))

    idx = { p:i for i, p in enumerate(profile.T[0])}
    cluster_idx = sorted([ [idx.get(c, -1), i] for i, c in enumerate(cluster.T[0]) if c in idx ])
    cluster = cluster[np.array(cluster_idx).T[1]]

    idx = { p:i for i, p in enumerate(cluster.T[0])}
    cluster_idx = sorted([ [idx.get(c, -1), i] for i, c in enumerate(profile.T[0]) if c in idx ])
    profile = profile[np.array(cluster_idx).T[1]]

    cluster.T[0] = np.arange(cluster.shape[0])
    profile.T[0] = np.arange(profile.shape[0])
    cluster = cluster.astype(int)
    profile = profile.astype(int)
    cluster[cluster < 0] = 0
    profile[profile < 0] = 0

    cluster = cluster[:, 1::stepwise]

    silhouette = get_silhouette(profile, cluster, stepwise, pool)
    similarity = get_similarity(normalized_mutual_info_score, cluster, stepwise, pool)

    with open(output+'.tsv', 'w') as fout:
        levels = ['HC{0}'.format(lvl*stepwise) for lvl in np.arange(silhouette.shape[0])]
        for lvl, ss in zip(levels, silhouette) :
            fout.write('#Silhouette\t{0}\t{1}\n'.format(lvl, ss))

        fout.write('\n#NMI\t{0}\n'.format('\t'.join(levels)))
        for lvl, nmis in zip(levels, similarity):
            fout.write('{0}\t{1}\n'.format(lvl, '\t'.join([ '{0:.3f}'.format(nmi) for nmi in nmis ])))
    fig, axs = plt.subplots(2, 2, \
                            figsize=(8, 12), \
                            gridspec_kw={'width_ratios':(12, 1),
                                         'height_ratios': (65, 35)})

    heatplot = axs[0, 0].imshow( (10*(np.log10(1-similarity))), \
                                norm=colors.TwoSlopeNorm(vmin=-30., vcenter=-10., vmax=0), \
                                cmap = 'RdBu',\
                                extent=[0, silhouette.shape[0]*stepwise, \
                                        silhouette.shape[0]*stepwise, 0])
    cb = fig.colorbar(heatplot, cax=axs[0, 1])
    axs[1, 0].plot(np.arange(silhouette.shape[0])*stepwise, silhouette,)
    axs[1, 0].set_xlim([0, silhouette.shape[0]*stepwise])
    axs[1, 1].remove()
    axs[0, 0].set_ylabel('HCs (allelic distances)')
    axs[0, 0].set_xlabel('HCs (allelic distances)')
    axs[1, 0].set_ylabel('Silhouette scores')
    axs[1, 0].set_xlabel('HCs (allelic distances)')
    cb.set_label('Normalized Mutual Information')
    cb.set_ticks([-30, -23.01, -20, -13.01, -10, -3.01, 0])
    cb.ax.set_yticklabels(['>=.999', '.995', '.99', '.95', '.9', '.5', '.0'])
    plt.savefig(output+'.pdf')
    pool.close()
    mapped_s = list(zip(levels, silhouette))
    mapped_sim = np.column_stack((names, cluster))
    return mapped_s, mapped_sim


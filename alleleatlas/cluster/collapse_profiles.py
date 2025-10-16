"""
Collapse similar samples in cgMLST profiles based on distance threshold.

This module provides functionality to:
1. Automatically determine a collapse threshold to reach a target number of clusters
2. Collapse samples with inter-sample distances below threshold
3. Preserve sample count information for downstream visualization (UMAP, MST, MDS)

Usage:
    from alleleatlas.cluster.collapse_profiles import auto_collapse_and_save
    
    collapsed_profile, metadata = auto_collapse_and_save(
        profile_path,
        outdir,
        target_nodes=10000,
        nproc=4
    )
    # Returns:
    #   collapsed_profile: path to collapsed .tsv with ST column as cluster ID
    #   metadata: dict with counts, thresholds, sample mappings
"""

from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.sparse.csgraph import minimum_spanning_tree

from alleleatlas.cluster.pHierCC import prepare_mat
from alleleatlas.cluster.getDistance import getDistance


def _find_collapse_threshold(mst_edges, n_samples, target_nodes):
    """
    Binary search to find collapse threshold that yields target number of clusters.
    
    Parameters:
        mst_edges: list of (u, v, weight) tuples from MST
        n_samples: total number of samples
        target_nodes: target number of clusters
    
    Returns:
        threshold: MST edge weight threshold that yields closest to target_nodes clusters
    """
    # Get sorted unique edge weights from MST
    weights = sorted(set(w for u, v, w in mst_edges))
    
    if not weights or len(weights) == 0:
        return 0
    
    # If target already less than samples, threshold must be high enough
    if target_nodes >= n_samples:
        return weights[-1] + 1  # Use a high threshold (no collapsing)
    
    # Binary search on weight values to find threshold that gives target_nodes
    best_threshold = weights[-1]
    best_diff = abs(n_samples - target_nodes)
    
    for threshold in weights:
        n_clusters = _count_clusters_at_threshold(mst_edges, n_samples, threshold)
        diff = abs(n_clusters - target_nodes)
        
        if diff < best_diff:
            best_diff = diff
            best_threshold = threshold
        
        # If we've crossed the target, check if we're getting worse
        if n_clusters <= target_nodes:
            break
    
    return best_threshold


def _count_clusters_at_threshold(mst_edges, n_samples, threshold):
    """Count number of clusters when collapsing edges < threshold."""
    parent = list(range(n_samples))
    
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    
    # Collapse edges below threshold
    for u, v, w in mst_edges:
        if w < threshold:
            union(u, v)
    
    # Count unique roots
    return len(set(find(i) for i in range(n_samples)))


def _collapse_samples(mst_edges, n_samples, threshold):
    """
    Collapse samples where MST edge < threshold.
    
    Returns:
        clusters: list of lists, each containing original sample indices
        cluster_map: dict mapping original sample index to cluster ID
    """
    parent = list(range(n_samples))
    
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    
    # Collapse edges below threshold
    for u, v, w in mst_edges:
        if w < threshold:
            union(u, v)
    
    # Group samples by cluster
    clusters_dict = defaultdict(list)
    for i in range(n_samples):
        root = find(i)
        clusters_dict[root].append(i)
    
    # Convert to list format and create mapping
    clusters = list(clusters_dict.values())
    cluster_map = {}
    for cluster_id, sample_indices in enumerate(clusters):
        for sample_idx in sample_indices:
            cluster_map[sample_idx] = cluster_id
    
    return clusters, cluster_map


def auto_collapse_and_save(profile_path, outdir, target_nodes=10000, nproc=4):
    """
    Automatically collapse a profile to target number of clusters and save.
    
    Parameters:
        profile_path (str): Path to input profile (.tsv, .csv, .gz, etc.)
        outdir (str): Output directory
        target_nodes (int): Target number of clusters (default 10000)
        nproc (int): Number of processes for distance computation
    
    Returns:
        tuple: (collapsed_profile_path, metadata_dict)
        
    metadata_dict contains:
        - n_original: original number of samples
        - n_collapsed: number of clusters
        - collapse_threshold: edge weight threshold used
        - sample_counts: dict mapping cluster ID to number of original samples
        - cluster_to_samples: dict mapping cluster ID to list of original sample indices
        - profile_path: path to collapsed profile
    """
    from multiprocessing import Pool
    
    profile_path = Path(profile_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Read and prepare matrix
    print(f"Loading profile from {profile_path}...")
    mat, names = prepare_mat(str(profile_path))
    n_samples = mat.shape[0]
    print(f"  {n_samples} samples, {mat.shape[1]-1} loci")
    
    # Check if collapse is needed
    if n_samples <= target_nodes:
        print(f"Profile has {n_samples} samples (< target {target_nodes}). No collapsing needed.")
        return str(profile_path), {
            'n_original': n_samples,
            'n_collapsed': n_samples,
            'collapse_threshold': 0,
            'sample_counts': {i: 1 for i in range(n_samples)},
            'cluster_to_samples': {i: [i] for i in range(n_samples)},
            'profile_path': str(profile_path),
        }
    
    # Compute distances
    print(f"Computing pairwise distances for {n_samples} samples...")
    pool = Pool(nproc)
    try:
        dist3 = getDistance(mat, 'dual_dist', pool, start=0, allowed_missing=0.03)
    finally:
        pool.close()
        pool.join()
    
    D = dist3[:, :, 0].astype(float)
    D = D + D.T
    np.fill_diagonal(D, 0.0)
    
    # Build MST
    print("Building minimum spanning tree...")
    mst_csr = minimum_spanning_tree(D, overwrite=False)
    mst_coo = mst_csr.tocoo()
    mst_edges = list(zip(mst_coo.row, mst_coo.col, mst_coo.data))
    
    # Find collapse threshold
    print(f"Finding collapse threshold to reach ~{target_nodes} clusters...")
    collapse_threshold = _find_collapse_threshold(mst_edges, n_samples, target_nodes)
    n_clusters = _count_clusters_at_threshold(mst_edges, n_samples, collapse_threshold)
    print(f"  Threshold: {collapse_threshold:.1f} -> {n_clusters} clusters")
    
    # Collapse samples
    print("Collapsing samples...")
    clusters, cluster_map = _collapse_samples(mst_edges, n_samples, collapse_threshold)
    
    # Build collapsed profile: use cluster ID as ST, merge loci by taking first sample
    print("Building collapsed profile...")
    collapsed_rows = []
    sample_counts = {}
    
    for cluster_id, sample_indices in enumerate(clusters):
        # Get representative (first sample in cluster)
        rep_idx = sample_indices[0]
        
        # Copy locus data from representative
        locus_data = mat[rep_idx, 1:].tolist()
        
        # Create row: [cluster_id, loci...]
        row = [str(cluster_id)] + [str(x) for x in locus_data]
        collapsed_rows.append(row)
        
        # Track count of original samples in this cluster
        sample_counts[cluster_id] = len(sample_indices)
    
    # Create DataFrame with same structure as input
    locus_cols = [f'L{i+1}' for i in range(mat.shape[1] - 1)]
    collapsed_df = pd.DataFrame(
        collapsed_rows,
        columns=['ST'] + locus_cols
    )
    
    # Save collapsed profile
    collapsed_profile_path = outdir / f"collapsed_{target_nodes}nodes_{collapse_threshold:.0f}thresh.tsv"
    collapsed_df.to_csv(collapsed_profile_path, sep='\t', index=False)
    print(f"Saved collapsed profile to: {collapsed_profile_path}")
    
    # Save metadata
    metadata = {
        'n_original': n_samples,
        'n_collapsed': n_clusters,
        'collapse_threshold': collapse_threshold,
        'target_nodes': target_nodes,
        'sample_counts': sample_counts,  # cluster_id -> count
        'cluster_to_samples': {i: cluster for i, cluster in enumerate(clusters)},  # cluster_id -> [sample indices]
        'profile_path': str(collapsed_profile_path),
    }
    
    # Save metadata as JSON
    import json
    metadata_path = outdir / f"collapse_metadata_{target_nodes}nodes_{collapse_threshold:.0f}thresh.json"
    # Convert to JSON-serializable format
    metadata_json = {
        'n_original': metadata['n_original'],
        'n_collapsed': metadata['n_collapsed'],
        'collapse_threshold': float(metadata['collapse_threshold']),
        'target_nodes': metadata['target_nodes'],
        'sample_counts': {str(k): v for k, v in metadata['sample_counts'].items()},
        'profile_path': metadata['profile_path'],
    }
    with open(metadata_path, 'w') as f:
        json.dump(metadata_json, f, indent=2)
    print(f"Saved metadata to: {metadata_path}")
    
    return str(collapsed_profile_path), metadata


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python -m alleleatlas.cluster.collapse_profiles <profile> [outdir] [target_nodes]")
        sys.exit(1)
    
    profile = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) > 2 else 'output_collapsed'
    target_nodes = int(sys.argv[3]) if len(sys.argv) > 3 else 10000
    
    collapsed_path, metadata = auto_collapse_and_save(profile, outdir, target_nodes, nproc=4)
    print(f"\n✓ Collapsed profile ready: {collapsed_path}")
    print(f"  Original samples: {metadata['n_original']}")
    print(f"  Collapsed clusters: {metadata['n_collapsed']}")
    print(f"  Collapse threshold: {metadata['collapse_threshold']:.1f}")

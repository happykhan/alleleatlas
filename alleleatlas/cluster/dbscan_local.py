"""DBSCAN-based local clustering for identifying tightly related samples.

This module provides local clustering capabilities to identify groups of
samples that are very similar (within certain distance thresholds), which
can be collapsed for computational efficiency.
"""

import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score, davies_bouldin_score
from rich.console import Console

console = Console()


def find_optimal_eps(distance_matrix, min_samples=5, method='silhouette', n_candidates=20):
    """Find optimal DBSCAN eps parameter via silhouette or k-distance method.
    
    Parameters:
        distance_matrix: Symmetric distance matrix (n_samples x n_samples) or 3D array from getDistance
        min_samples: Minimum samples per cluster
        method: 'silhouette' (default) or 'kdistance'
        n_candidates: Number of eps values to test
        
    Returns:
        (optimal_eps, scores_dict)
    """
    # Handle 3D distance matrix from getDistance (extract first channel - the distances)
    # Shape is [n-start, n, 2] representing lower triangle of distance matrix
    if distance_matrix.ndim == 3:
        dist_2d = distance_matrix[:, :, 0].astype(float)
        # Reconstruct symmetric matrix
        n_samples = dist_2d.shape[1]
        start = n_samples - dist_2d.shape[0]
        symmetric_dist = np.zeros((n_samples, n_samples), dtype=float)
        symmetric_dist[start:, :] = dist_2d
        # Mirror to upper triangle
        for i in range(n_samples):
            for j in range(i):
                symmetric_dist[j, i] = symmetric_dist[i, j]
        distance_matrix = symmetric_dist
    
    # Convert to condensed distance for DBSCAN
    from scipy.spatial.distance import squareform
    distances_condensed = squareform(distance_matrix)
    
    if method == 'kdistance':
        # K-distance method: find "knee" in sorted k-distances
        k_distances = np.array([distance_matrix[i].mean() for i in range(len(distance_matrix))])
        k_distances = np.sort(k_distances)
        
        # Simple knee detection: find biggest gap
        diffs = np.diff(k_distances)
        knee_idx = np.argmax(diffs)
        optimal_eps = k_distances[knee_idx]
        
        return optimal_eps, {'method': 'kdistance', 'knee_idx': knee_idx}
    
    else:  # silhouette method
        console.print("  [dim]Testing DBSCAN eps values via silhouette...[/dim]")
        
        # Test range of eps values
        eps_range = np.linspace(
            np.percentile(distances_condensed, 10),
            np.percentile(distances_condensed, 90),
            n_candidates
        )
        
        best_eps = None
        best_score = -1
        scores = {}
        
        for eps in eps_range:
            try:
                clustering = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
                labels = clustering.fit_predict(distance_matrix)
                
                n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
                n_noise = list(labels).count(-1)
                
                # Skip if too few clusters or too much noise
                if n_clusters < 2 or n_noise > 0.5 * len(labels):
                    continue
                
                sil_score = silhouette_score(distance_matrix, labels, metric='precomputed')
                scores[eps] = {
                    'silhouette': sil_score,
                    'n_clusters': n_clusters,
                    'n_noise': n_noise
                }
                
                if sil_score > best_score:
                    best_score = sil_score
                    best_eps = eps
                    
            except Exception as e:
                console.print(f"    [dim]eps={eps:.2f}: failed ({type(e).__name__})[/dim]", overflow="fold")
                continue
        
        if best_eps is None:
            # Fallback to median distance
            best_eps = np.median(distances_condensed)
            console.print(f"  [yellow]⚠ No good eps found, using median distance: {best_eps:.2f}[/yellow]")
        else:
            console.print(f"  [green]✓[/green] Optimal eps: {best_eps:.2f} (silhouette: {best_score:.3f})")
        
        return best_eps, scores


def run_dbscan_clustering(distance_matrix, eps=None, min_samples=5, auto_eps=True):
    """Run DBSCAN clustering on distance matrix.
    
    Parameters:
        distance_matrix: Symmetric distance matrix (n_samples x n_samples) or 3D array from getDistance
        eps: Distance threshold (if None and auto_eps=True, find optimal)
        min_samples: Minimum samples per cluster
        auto_eps: If True, find optimal eps when eps=None
        
    Returns:
        (labels, eps_used, metrics_dict)
    """
    # Handle 3D distance matrix from getDistance (extract first channel - the distances)
    # Shape is [n-start, n, 2] representing lower triangle of distance matrix
    if distance_matrix.ndim == 3:
        dist_2d = distance_matrix[:, :, 0].astype(float)
        # Reconstruct symmetric matrix
        n_samples = dist_2d.shape[1]
        start = n_samples - dist_2d.shape[0]
        symmetric_dist = np.zeros((n_samples, n_samples), dtype=float)
        symmetric_dist[start:, :] = dist_2d
        # Mirror to upper triangle
        for i in range(n_samples):
            for j in range(i):
                symmetric_dist[j, i] = symmetric_dist[i, j]
        distance_matrix = symmetric_dist
    
    if eps is None and auto_eps:
        eps, _ = find_optimal_eps(distance_matrix, min_samples=min_samples, method='silhouette')
    elif eps is None:
        eps = 10  # reasonable default
    
    console.print(f"Running DBSCAN with eps={eps:.2f}, min_samples={min_samples}")
    
    # Run clustering
    clustering = DBSCAN(eps=eps, min_samples=min_samples, metric='precomputed')
    labels = clustering.fit_predict(distance_matrix)
    
    # Compute metrics
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = list(labels).count(-1)
    n_samples = len(labels)
    
    metrics = {
        'eps': eps,
        'n_clusters': n_clusters,
        'n_noise': n_noise,
        'n_samples': n_samples,
        'noise_fraction': n_noise / n_samples if n_samples > 0 else 0,
    }
    
    # Silhouette score (only for non-noise points if there are any)
    valid_mask = labels != -1
    if valid_mask.sum() > 1 and n_clusters > 1:
        try:
            # Convert to Euclidean distance for silhouette score
            sil_score = silhouette_score(
                distance_matrix,
                labels,
            )
            metrics['silhouette_score'] = sil_score
        except Exception:
            pass
    
    # Davies-Bouldin index (lower is better)
    if n_clusters > 1 and valid_mask.sum() > 0:
        try:
            db_score = davies_bouldin_score(
                distance_matrix,
                labels,
            )
            metrics['davies_bouldin_index'] = db_score
        except Exception:
            pass
    
    console.print(f"  [green]✓[/green] Found {n_clusters} clusters with {n_noise} noise points")
    if 'silhouette_score' in metrics:
        console.print(f"  Silhouette score: {metrics['silhouette_score']:.3f}")
    
    return labels, eps, metrics


def get_cluster_sizes(labels):
    """Get size of each cluster (excluding noise).
    
    Parameters:
        labels: DBSCAN cluster labels (including -1 for noise)
        
    Returns:
        dict: {cluster_id: count}
    """
    unique, counts = np.unique(labels, return_counts=True)
    size_dict = {}
    for uid, count in zip(unique, counts):
        if uid != -1:  # Skip noise
            size_dict[uid] = count
    return size_dict


__all__ = [
    'find_optimal_eps',
    'run_dbscan_clustering',
    'get_cluster_sizes',
]

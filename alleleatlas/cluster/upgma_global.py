"""UPGMA-based global clustering for identifying broader structural groups.

This module provides hierarchical clustering using UPGMA (Unweighted Pair Group
Method with Arithmetic Mean) to discover the overall structure and generate
dendrograms for visualization.
"""

import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import matplotlib.pyplot as plt
from pathlib import Path
from rich.console import Console

console = Console()


def compute_upgma_linkage(distance_matrix, method='single'):
    """Compute UPGMA hierarchical clustering linkage.
    
    Parameters:
        distance_matrix: Symmetric distance matrix (n_samples x n_samples) or 3D array from getDistance
        method: Linkage method ('single', 'complete', 'average', 'weighted', 'centroid', etc.)
                Default 'single' for single-linkage clustering
        
    Returns:
        Z: Linkage matrix from scipy.cluster.hierarchy.linkage
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
    
    console.print(f"Computing {method}-linkage hierarchical clustering...")
    
    # Convert to condensed distance form
    distances_condensed = squareform(distance_matrix)
    
    # Compute linkage
    Z = linkage(distances_condensed, method=method)
    
    console.print(f"  [green]✓[/green] Linkage computed for {distance_matrix.shape[0]} samples")
    
    return Z


def get_clusters_at_distance(Z, distance_threshold):
    """Get cluster labels by cutting dendrogram at specific distance threshold.
    
    Parameters:
        Z: Linkage matrix
        distance_threshold: Distance at which to cut the dendrogram
        
    Returns:
        labels: Cluster labels for each sample
    """
    labels = fcluster(Z, distance_threshold, criterion='distance')
    return labels


def get_clusters_at_level(Z, n_clusters):
    """Get cluster labels by specifying number of clusters.
    
    Parameters:
        Z: Linkage matrix
        n_clusters: Number of clusters to extract
        
    Returns:
        labels: Cluster labels for each sample
    """
    labels = fcluster(Z, n_clusters, criterion='maxclust')
    return labels


def plot_dendrogram(Z, sample_names=None, max_d=None, figsize=(16, 8), outpath=None):
    """Plot and save dendrogram from linkage matrix.
    
    Parameters:
        Z: Linkage matrix
        sample_names: Optional list of sample names (if None, use indices)
        max_d: Optional vertical line to draw at this distance
        figsize: Figure size (width, height)
        outpath: Path to save figure (if None, don't save)
        
    Returns:
        dendro: Dendrogram dict from scipy
    """
    console.print("Plotting dendrogram...")
    
    plt.figure(figsize=figsize)
    
    kwargs = {
        'leaf_rotation': 90,
        'leaf_font_size': 8,
        'no_labels': sample_names is None and len(Z) > 100,  # Hide labels if too many
    }
    
    dendro = dendrogram(Z, **kwargs)
    
    plt.xlabel('Sample Index' if sample_names is None else 'Sample')
    plt.ylabel('Distance')
    plt.title(f'UPGMA Hierarchical Clustering Dendrogram ({Z.shape[0] + 1} samples)')
    
    if max_d is not None:
        plt.axhline(y=max_d, c='red', linestyle='--', linewidth=2, label=f'Cut at d={max_d:.2f}')
        plt.legend()
    
    if outpath is not None:
        Path(outpath).parent.mkdir(parents=True, exist_ok=True)
        try:
            # Try to save with minimal formatting
            plt.savefig(outpath, dpi=100, pad_inches=0.2)
            console.print(f"  [green]✓[/green] Dendrogram saved to {outpath}")
        except Exception as e:
            console.print(f"  [yellow]⚠[/yellow] Could not save dendrogram: {e}")
    
    plt.close()
    return dendro


def get_optimal_distance_levels(Z, n_levels=5):
    """Get suggested distance levels from dendrogram for clustering.
    
    Finds the n largest jumps in the dendrogram's distance values, which
    typically correspond to meaningful cluster merges.
    
    Parameters:
        Z: Linkage matrix
        n_levels: Number of suggested levels to return
        
    Returns:
        list: Suggested distance thresholds (sorted ascending)
    """
    # Last column of Z contains distances
    distances = Z[:, 2]
    
    # Find the largest jumps (differences)
    diffs = np.diff(distances)
    largest_jumps = np.argsort(diffs)[-n_levels:][::-1]
    
    # Get the distances at those jumps
    suggested_distances = distances[largest_jumps]
    suggested_distances = np.sort(suggested_distances)
    
    return suggested_distances.tolist()


def analyze_dendrogram(Z, distance_matrix):
    """Analyze dendrogram structure and return statistics.
    
    Parameters:
        Z: Linkage matrix
        distance_matrix: Original distance matrix
        
    Returns:
        dict: Analysis results including:
            - n_samples: Number of samples
            - n_merges: Number of merge steps
            - min_distance: Minimum distance in linkage
            - max_distance: Maximum distance in linkage
            - mean_distance: Mean distance
            - suggested_levels: List of suggested cut distances
    """
    distances = Z[:, 2]
    
    analysis = {
        'n_samples': distance_matrix.shape[0],
        'n_merges': len(Z),
        'min_distance': float(distances.min()),
        'max_distance': float(distances.max()),
        'mean_distance': float(distances.mean()),
        'median_distance': float(np.median(distances)),
        'suggested_levels': get_optimal_distance_levels(Z, n_levels=5),
    }
    
    return analysis


def get_clusters_at_multiple_levels(Z, distance_levels):
    """Get cluster assignments at multiple distance levels.
    
    Parameters:
        Z: Linkage matrix
        distance_levels: List of distance thresholds
        
    Returns:
        dict: {distance: labels_array}
    """
    results = {}
    for dist in distance_levels:
        labels = get_clusters_at_distance(Z, dist)
        results[dist] = labels
    return results


def run_upgma_clustering(distance_matrix, sample_names=None, outdir=None):
    """Complete UPGMA clustering workflow.
    
    Parameters:
        distance_matrix: Symmetric distance matrix
        sample_names: Optional sample names/IDs
        outdir: Output directory for plots
        
    Returns:
        dict: Results containing:
            - Z: Linkage matrix
            - analysis: Dendrogram analysis
            - dendro_plot_path: Path to saved dendrogram
            - suggested_levels: Suggested distance cutoffs
    """
    # Compute linkage
    Z = compute_upgma_linkage(distance_matrix)
    
    # Analyze structure
    analysis = analyze_dendrogram(Z, distance_matrix)
    
    # Plot dendrogram
    dendro_plot_path = None
    if outdir is not None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        dendro_plot_path = str(outdir / 'upgma_dendrogram.png')
        
        # Determine reasonable max_d for visualization
        levels = analysis['suggested_levels']
        max_d_viz = levels[-1] if levels else None
        
        plot_dendrogram(
            Z,
            sample_names=sample_names,
            max_d=max_d_viz,
            outpath=dendro_plot_path
        )
    
    # Print analysis
    console.print("\n[bold]UPGMA Analysis:[/bold]")
    console.print(f"  Samples: {analysis['n_samples']}")
    console.print(f"  Distance range: {analysis['min_distance']:.2f} - {analysis['max_distance']:.2f}")
    console.print(f"  Median distance: {analysis['median_distance']:.2f}")
    console.print(f"  Suggested levels: {[f'{d:.2f}' for d in analysis['suggested_levels'][:3]]}")
    
    return {
        'Z': Z,
        'analysis': analysis,
        'dendro_plot_path': dendro_plot_path,
        'suggested_levels': analysis['suggested_levels'],
    }


__all__ = [
    'compute_upgma_linkage',
    'get_clusters_at_distance',
    'get_clusters_at_level',
    'plot_dendrogram',
    'get_optimal_distance_levels',
    'analyze_dendrogram',
    'get_clusters_at_multiple_levels',
    'run_upgma_clustering',
]

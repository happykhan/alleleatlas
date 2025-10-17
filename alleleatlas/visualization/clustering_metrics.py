"""Clustering validation and visualization functions.

Provides MDS visualization, silhouette scoring, and clustering metrics plots
for validating DBSCAN clustering results.
"""

from pathlib import Path
import numpy as np
import matplotlib
import os 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from sklearn.metrics import silhouette_score
from rich.console import Console

console = Console()


def _symmetrize_distance_matrix_from_phiercc(dist_3d):
    """Convert 3D lower-triangle distance matrix (from pHierCC getDistance) to symmetric 2D.
    
    pHierCC's getDistance returns dist with shape (n_samples-start, n_samples, 2) where:
    - dist[i-start, j, :] contains the distance between sample i and sample j
    - Only the lower triangle is populated (j < i)
    - dist[:, :, 0] is the "dual" distance (with missing data correction)
    - dist[:, :, 1] is the "query" distance (query-specific missing penalty)
    
    This function reconstructs the full symmetric matrix by:
    1. Taking only the dual distances [:, :, 0]
    2. Properly placing them in the full matrix accounting for the row offset
    3. Mirroring across the diagonal to create a symmetric matrix
    
    Parameters:
        dist_3d: 3D array of shape (n_samples-start, n_samples, 2) from getDistance
        
    Returns:
        2D symmetric distance matrix of shape (n_samples, n_samples)
    """
    if dist_3d.ndim != 3 or dist_3d.shape[2] != 2:
        raise ValueError(f"Expected 3D array with shape (n, m, 2), got {dist_3d.shape}")
    
    n_query, n_samples, _ = dist_3d.shape
    start = n_samples - n_query  # Offset for the first query sample
    
    # Extract dual distances (the main clustering-relevant distance metric)
    dist_dual = dist_3d[:, :, 0].astype(float)
    
    # Create full symmetric matrix
    symmetric = np.zeros((n_samples, n_samples), dtype=float)
    
    # Place the lower-triangle data into the matrix
    for i in range(start, n_samples):
        for j in range(i):
            symmetric[i, j] = dist_dual[i - start, j]
    
    # Mirror to create symmetric matrix
    for i in range(n_samples):
        for j in range(i):
            symmetric[j, i] = symmetric[i, j]
    
    return symmetric


def plot_mds_with_clusters(distance_matrix, labels, outdir):
    """Plot MDS embedding with DBSCAN cluster colors.
    
    Creates a 2D projection of the samples using Multidimensional Scaling,
    with points colored by their DBSCAN cluster assignment. Noise points
    (cluster -1) are shown in gray.
    
    Parameters:
        distance_matrix: Distance matrix (2D symmetric or 3D lower-triangle)
        labels: DBSCAN cluster labels (array-like)
        outdir: Output directory (str or Path)
        
    Returns:
        True if successful, False if error occurred
        
    Notes:
        - Automatically handles 3D lower-triangle format from getDistance()
        - Uses scikit-learn MDS with precomputed distances
        - Saves output to: outdir/mds_clusters.png
    """
    try:
        outdir = Path(outdir)
        
        # Handle 3D distance matrix (lower triangle format from pHierCC)
        if distance_matrix.ndim == 3:
            distance_matrix = _symmetrize_distance_matrix_from_phiercc(distance_matrix)
        
        # Compute MDS
        mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
        embedding = mds.fit_transform(distance_matrix)
        
        # Plot with cluster colors
        plt.figure(figsize=(10, 8))
        
        # Plot each cluster with different color
        unique_labels = np.unique(labels)
        colors = plt.cm.get_cmap('tab10')(np.linspace(0, 1, len(unique_labels)))
        
        for label, color in zip(unique_labels, colors):
            mask = labels == label
            if label == -1:
                # Noise points
                plt.scatter(embedding[mask, 0], embedding[mask, 1], 
                           c='gray', s=30, alpha=0.5, label='Noise')
            else:
                plt.scatter(embedding[mask, 0], embedding[mask, 1], 
                           c=[color], s=50, alpha=0.7, label=f'Cluster {int(label)}')
        
        plt.xlabel('MDS1')
        plt.ylabel('MDS2')
        plt.title('MDS Plot with DBSCAN Clusters')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        plt.tight_layout()
        plt.savefig(outdir / 'mds_clusters.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        console.print('  ✓ MDS plot saved')
        return True
        
    except Exception as exc:  # pylint: disable=broad-except
        console.print(f'  [yellow]⚠[/yellow] Could not create MDS plot: {exc}')
        return False


def plot_silhouette_scores(distance_matrix, labels, outdir):
    """Plot silhouette coefficient for clustering validation.
    
    Creates a bar chart showing the silhouette coefficient, which measures
    how well-separated the clusters are. Range [-1, 1] where:
    - 1.0 = perfect clustering
    - 0.5 = good clustering
    - 0.0 = random clustering
    - <0.0 = overlapping clusters
    
    Parameters:
        distance_matrix: Distance matrix (2D symmetric or 3D lower-triangle)
        labels: Cluster labels (array-like)
        outdir: Output directory (str or Path)
        
    Returns:
        True if successful, False if error occurred
        
    Notes:
        - Only plots if there are 2+ clusters
        - Uses precomputed distance metric
        - Methodology matches HCCeval.py approach
        - Saves output to: outdir/silhouette_score.png
    """
    try:
        outdir = Path(outdir)
        
        # Handle 3D distance matrix (lower triangle format from pHierCC)
        if distance_matrix.ndim == 3:
            distance_matrix = _symmetrize_distance_matrix_from_phiercc(distance_matrix)
        
        # Compute silhouette scores
        n_clusters = len(np.unique(labels)) - (1 if -1 in labels else 0)
        
        if n_clusters > 1:
            score = silhouette_score(distance_matrix, labels, metric='precomputed')
            
            plt.figure(figsize=(8, 6))
            plt.bar(range(1), [float(score)], color='steelblue', alpha=0.7)
            plt.ylabel('Score')
            plt.title(f'Silhouette Score (Clusters: {n_clusters})')
            plt.ylim([-1.0, 1.0])
            plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
            plt.text(0, float(score) + 0.05, f'{score:.3f}', ha='center', fontweight='bold')
            plt.xticks(range(1), ['Silhouette Score'])
            plt.tight_layout()
            plt.savefig(outdir / 'silhouette_score.png', dpi=150)
            plt.close()
            
            console.print(f'  ✓ Silhouette plot saved (score: {score:.3f})')
            return True
        else:
            console.print('  [yellow]⚠[/yellow] Need 2+ clusters to compute silhouette score')
            return False
            
    except Exception as exc:  # pylint: disable=broad-except
        console.print(f'  [yellow]⚠[/yellow] Could not create silhouette plot: {exc}')
        return False


def plot_clustering_metrics(distance_matrix, labels, outdir):
    """Plot clustering validation metrics in a multi-panel figure.
    
    Creates a 2-panel comparison plot showing:
    - Left: Silhouette coefficient (clustering quality metric)
    - Right: Cluster size distribution (cluster balance)
    
    This approach matches the methodology used in HCCeval.py for
    comprehensive clustering assessment.
    
    Parameters:
        distance_matrix: Distance matrix (2D symmetric or 3D lower-triangle)
        labels: Cluster labels (array-like)
        outdir: Output directory (str or Path)
        
    Returns:
        True if successful, False if error occurred
        
    Notes:
        - Automatically handles 3D lower-triangle format from getDistance()
        - Only computes silhouette if 2+ clusters present
        - Saves output to: outdir/clustering_metrics.png
        - Uses different colors for visual distinction
    """
    try:
        outdir = Path(outdir)
        
        # Handle 3D distance matrix (lower triangle format from pHierCC)
        if distance_matrix.ndim == 3:
            distance_matrix = _symmetrize_distance_matrix_from_phiercc(distance_matrix)
        
        n_clusters = len(np.unique(labels)) - (1 if -1 in labels else 0)
        
        # Compute metrics
        metrics = {}
        
        if n_clusters > 1:
            sil_score = silhouette_score(distance_matrix, labels, metric='precomputed')
            metrics['Silhouette'] = float(sil_score)
        
        # Create comparison plots
        axes_array = plt.subplots(1, 2, figsize=(12, 5))[1]
        
        # Plot 1: Silhouette score
        if 'Silhouette' in metrics:
            axes_array[0].bar(range(1), [metrics['Silhouette']], color='steelblue', alpha=0.7)
            axes_array[0].set_ylabel('Score')
            axes_array[0].set_title('Silhouette Score')
            axes_array[0].set_ylim([-1.0, 1.0])
            axes_array[0].axhline(y=0, color='black', linestyle='-', linewidth=0.5)
            axes_array[0].text(0, metrics['Silhouette'] + 0.05, 
                        f'{metrics["Silhouette"]:.3f}', ha='center', fontweight='bold')
            axes_array[0].set_xticks(range(1))
            axes_array[0].set_xticklabels(['Score'])
        
        # Plot 2: Cluster size distribution
        unique_labels, counts = np.unique(labels, return_counts=True)
        cluster_names = [f'Cluster {lab}' if lab != -1 else 'Noise' for lab in unique_labels]
        axes_array[1].bar(cluster_names, counts, color='coral', alpha=0.7)
        axes_array[1].set_ylabel('Number of Samples')
        axes_array[1].set_title('Cluster Size Distribution')
        axes_array[1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(outdir / 'clustering_metrics.png', dpi=150)
        plt.close()
        
        console.print('  ✓ Clustering metrics plot saved')
        return True
        
    except Exception as exc:  # pylint: disable=broad-except
        console.print(f'  [yellow]⚠[/yellow] Could not create metrics plot: {exc}')
        return False


def plot_hcceval_style(distance_matrix, labels, outdir, title="Clustering Evaluation"):
    """Plot HCCeval-style 2-panel figure with NMI heatmap and silhouette scores.
    
    Creates a 2-panel figure matching HCCeval.py output style:
    - Top panel: Heatmap of normalized mutual information (NMI) between clustering levels
    - Bottom panel: Silhouette scores across hierarchical levels
    
    This provides a comprehensive view of clustering stability and quality
    across multiple distance thresholds, matching the pHierCC evaluation methodology.
    
    Parameters:
        distance_matrix: Distance matrix (2D symmetric or 3D lower-triangle)
        labels: Cluster labels (array-like, typically from DBSCAN or hierarchical clustering)
        outdir: Output directory (str or Path)
        title: Title for the figure (optional)
        
    Returns:
        True if successful, False if error occurred
        
    Notes:
        - Automatically handles 3D lower-triangle format from getDistance()
        - Requires 2+ clusters for meaningful silhouette computation
        - Uses normalized mutual information (NMI) for level similarity
        - Saves output to: outdir/hcceval_evaluation.png
        - Color scale follows HCCeval conventions (log-transformed NMI)
    """
    try:
        outdir = Path(outdir)
        
        # Handle 3D distance matrix (lower triangle format from pHierCC)
        if distance_matrix.ndim == 3:
            distance_matrix = _symmetrize_distance_matrix_from_phiercc(distance_matrix)
        
        # Ensure it's a valid square symmetric matrix
        if distance_matrix.ndim != 2 or distance_matrix.shape[0] != distance_matrix.shape[1]:
            console.print(f'  [yellow]⚠[/yellow] Invalid distance matrix shape: {distance_matrix.shape}')
            return False
        
        n_clusters = len(np.unique(labels)) - (1 if -1 in labels else 0)
        
        if n_clusters < 2:
            console.print('  [yellow]⚠[/yellow] Need 2+ clusters for HCCeval-style evaluation')
            return False
        
        # Compute silhouette score for the current clustering
        silhouette_score_val = silhouette_score(distance_matrix, labels, metric='precomputed')
        
        # Create figure with subplots matching HCCeval style
        # Using a 2x2 grid where we remove the bottom-right
        fig = plt.figure(figsize=(10, 10))
        gs = fig.add_gridspec(2, 2, width_ratios=(12, 1), height_ratios=(2, 1),
                              hspace=0.3, wspace=0.3)
        
        ax_heatmap = fig.add_subplot(gs[0, 0])
        ax_cbar = fig.add_subplot(gs[0, 1])
        ax_silhouette = fig.add_subplot(gs[1, 0])
        
        # ========== Top Panel: NMI Heatmap ==========
        # Create a simple NMI matrix showing label agreement at different "thresholds"
        # For now, we'll create a synthetic multi-level evaluation by sampling at different distance percentiles
        
        n_levels = min(10, len(np.unique(labels)) + 3)  # Create 3-10 evaluation levels
        
        # Generate multiple clusterings at different distance thresholds
        # This simulates hierarchical clustering at different levels
        from sklearn.cluster import AgglomerativeClustering
        from sklearn.metrics import normalized_mutual_info_score
        
        multi_labels = []
        distance_condensed = distance_matrix.copy()
        
        # Create clusterings at different levels
        for n_clust in range(2, min(n_levels + 2, len(labels))):
            try:
                hierarchical = AgglomerativeClustering(
                    n_clusters=n_clust, 
                    linkage='average',
                    metric='precomputed'
                )
                hier_labels = hierarchical.fit_predict(distance_condensed)
                multi_labels.append(hier_labels)
            except Exception:
                break
        
        if len(multi_labels) < 2:
            # Fallback: just use current labels repeated
            multi_labels = [labels] * min(3, len(labels))
        
        # Compute NMI matrix between all pairs of clustering levels
        nmi_matrix = np.ones((len(multi_labels), len(multi_labels)), dtype=np.float64)
        for i, labels_i in enumerate(multi_labels):
            for j, labels_j in enumerate(multi_labels[i+1:], start=i+1):
                try:
                    nmi_val = normalized_mutual_info_score(labels_i, labels_j)
                    nmi_matrix[i, j] = nmi_val
                    nmi_matrix[j, i] = nmi_val
                except Exception:
                    nmi_matrix[i, j] = 0.5
                    nmi_matrix[j, i] = 0.5
        
        # Clamp values to valid range for log transformation
        nmi_matrix_safe = np.clip(nmi_matrix, 0.001, 0.999)
        
        # Apply log transformation matching HCCeval's color scale
        heatmap_data = 10 * (np.log10(1 - nmi_matrix_safe))
        
        # Create heatmap
        im = ax_heatmap.imshow(
            heatmap_data,
            cmap='RdBu',
            extent=(0, len(multi_labels), len(multi_labels), 0)
        )
        
        ax_heatmap.set_ylabel('Clustering Levels')
        ax_heatmap.set_xlabel('Clustering Levels')
        ax_heatmap.set_title(f'{title} - Multi-Level Clustering Similarity')
        
        # Colorbar
        cbar = fig.colorbar(im, cax=ax_cbar, label='NMI (log scale)')
        cbar.set_ticks([-30, -23.01, -20, -13.01, -10, -3.01, 0])
        cbar.ax.set_yticklabels(['≥.999', '.995', '.99', '.95', '.9', '.5', '.0'])
        
        # ========== Bottom Panel: Silhouette Scores ==========
        # Compute silhouette scores for each clustering level
        silhouette_scores_arr = []
        for hier_labels in multi_labels:
            try:
                sil_score = silhouette_score(distance_matrix, hier_labels, metric='precomputed')
                silhouette_scores_arr.append(float(sil_score))
            except Exception:
                silhouette_scores_arr.append(0.0)
        
        x_levels = np.arange(len(silhouette_scores_arr))
        ax_silhouette.plot(x_levels, silhouette_scores_arr, 'o-', linewidth=2, markersize=6, color='steelblue')
        ax_silhouette.set_xlim((-0.5, len(silhouette_scores_arr) - 0.5))
        ax_silhouette.set_ylim((-1.0, 1.0))
        ax_silhouette.set_ylabel('Silhouette Score')
        ax_silhouette.set_xlabel('Clustering Level')
        ax_silhouette.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax_silhouette.grid(True, alpha=0.3)
        
        # Add current score annotation
        ax_silhouette.text(0.5, 0.95, f'Current clustering silhouette: {silhouette_score_val:.3f}',
                          transform=ax_silhouette.transAxes,
                          ha='center', va='top',
                          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.savefig(outdir / 'hcceval_evaluation.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        console.print('  ✓ HCCeval-style evaluation plot saved')
        return True
        
    except Exception as exc:  # pylint: disable=broad-except
        console.print(f'  [yellow]⚠[/yellow] Could not create HCCeval plot: {exc}')
        import traceback
        traceback.print_exc()
        return False


def plot_hiercc_group_sizes(hiercc_groups, stepwise, output_dir):
    """Plot the number of groups at each hierarchical level.
    
    Parameters:
        hiercc_groups: 2D array of cluster assignments (n_samples × n_levels)
        stepwise: Step size between levels
        output_dir: Directory to save plot
    """
    import matplotlib.pyplot as plt
    
    # Count unique groups at each level
    group_counts = []
    distances = []
    
    for level_idx, groups in enumerate(hiercc_groups.T):
        n_groups = len(set(groups[groups > 0]))  # Count unique non-zero groups
        distance = level_idx * stepwise
        group_counts.append(n_groups)
        distances.append(distance)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(distances, group_counts, marker='o', linewidth=2, markersize=4, color='steelblue')
    ax.set_xlabel('Hierarchical Distance (allelic differences)', fontsize=12)
    ax.set_ylabel('Number of Clusters', fontsize=12)
    ax.set_title('HierCC: Number of Clusters at Each Level', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Add threshold line if there's a natural collapse point
    collapse_idx = None
    for idx, count in enumerate(group_counts):
        if count < 1000:
            collapse_idx = idx
            break
    
    if collapse_idx is not None and collapse_idx > 0 :
        ax.axvline(distances[collapse_idx], color='red', linestyle='--', alpha=0.7, 
                    label=f'Suggested threshold: {distances[collapse_idx]}')
        ax.legend()
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, 'hiercc_group_counts.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    console.print(f"  ✓ Saved group count plot to {output_path}")
    return output_path, group_counts

__all__ = [
    'plot_mds_with_clusters',
    'plot_silhouette_scores',
    'plot_clustering_metrics',
    'plot_hcceval_style',
]

"""Permutation testing for statistical validation of clustering results.

Generates null distributions by permuting cluster labels and comparing
real clustering metrics to randomly shuffled versions to assess significance.
"""

import numpy as np
from sklearn.metrics import silhouette_score
from multiprocessing import Pool
from rich.console import Console
from rich.progress import Progress

console = Console()


def compute_clustering_metrics(distance_matrix, labels, metric_types=None):
    """Compute clustering quality metrics.
    
    Parameters:
        distance_matrix: Symmetric distance matrix or 3D array from getDistance
        labels: Cluster labels
        metric_types: List of metrics to compute (default: ['silhouette', 'nmi_self'])
        
    Returns:
        dict: Computed metrics
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
    
    if metric_types is None:
        metric_types = ['silhouette']
    
    metrics = {}
    
    if 'silhouette' in metric_types and len(np.unique(labels)) > 1:
        try:
            sil = silhouette_score(distance_matrix, labels, metric='precomputed')
            metrics['silhouette'] = sil
        except Exception as e:
            console.print(f"  [yellow]Warning: silhouette computation failed: {e}[/yellow]")
    
    if 'davies_bouldin' in metric_types and len(np.unique(labels)) > 1:
        try:
            from sklearn.metrics import davies_bouldin_score
            db = davies_bouldin_score(distance_matrix, labels)
            metrics['davies_bouldin'] = db
        except Exception:
            pass
    
    return metrics


def permute_labels(labels, n_permutations=1, random_state=None):
    """Generate random permutations of labels.
    
    Parameters:
        labels: Original cluster labels
        n_permutations: Number of permuted versions to generate
        random_state: Seed for reproducibility
        
    Returns:
        list: List of permuted label arrays
    """
    if random_state is not None:
        np.random.seed(random_state)
    
    permuted = []
    for _ in range(n_permutations):
        perm_labels = np.random.permutation(labels)
        permuted.append(perm_labels)
    
    return permuted


def _compute_metrics_on_permutation(args):
    """Worker function for parallel permutation testing.
    
    Args:
        (distance_matrix, perm_labels, metric_types)
    
    Returns:
        dict of metric values
    """
    distance_matrix, perm_labels, metric_types = args
    try:
        return compute_clustering_metrics(distance_matrix, perm_labels, metric_types)
    except Exception as e:
        console.print(f"[yellow]Warning: permutation failed: {e}[/yellow]")
        return {}


def run_permutation_tests(
    distance_matrix,
    real_labels,
    n_permutations=100,
    metric_types=None,
    n_procs=4,
    random_state=None
):
    """Run permutation testing on clustering results.
    
    Parameters:
        distance_matrix: Symmetric distance matrix
        real_labels: Real cluster labels from clustering
        n_permutations: Number of permutations to test (default: 100)
        metric_types: Which metrics to compute (default: ['silhouette'])
        n_procs: Number of parallel processes
        random_state: Seed for reproducibility
        
    Returns:
        dict: Results including:
            - real_metrics: Metrics from real clustering
            - null_distributions: Dict of permutation metric distributions
            - p_values: Dict of p-values (fraction of permutations >= real value)
            - effect_sizes: Effect sizes (Cohen's d)
    """
    if metric_types is None:
        metric_types = ['silhouette']
    
    console.print("\n[bold]Permutation Testing:[/bold]")
    console.print(f"  Real labels: {len(np.unique(real_labels))} clusters")
    console.print(f"  Permutations: {n_permutations}")
    console.print(f"  Processes: {n_procs}")
    
    # Compute metrics on real clustering
    real_metrics = compute_clustering_metrics(distance_matrix, real_labels, metric_types)
    # Print metrics nicely (extract float values)
    metrics_str = ", ".join([f"{k}: {float(v):.4f}" for k, v in real_metrics.items()])
    console.print(f"  Real metrics: {metrics_str}")
    
    # Generate permutations
    if random_state is not None:
        np.random.seed(random_state)
    
    permuted_labels = permute_labels(real_labels, n_permutations, random_state=random_state)
    
    # Prepare work items
    work_items = [
        (distance_matrix, perm_labels, metric_types)
        for perm_labels in permuted_labels
    ]
    
    # Run in parallel
    null_distributions = {m: [] for m in metric_types}
    
    with Progress() as progress:
        task = progress.add_task("[cyan]Computing permutation metrics...", total=n_permutations)
        
        with Pool(n_procs) as pool:
            for perm_metrics in pool.imap_unordered(_compute_metrics_on_permutation, work_items):
                for metric_name, value in perm_metrics.items():
                    null_distributions[metric_name].append(value)
                progress.update(task, advance=1)
    
    # Convert lists to arrays
    for metric_name in null_distributions:
        null_distributions[metric_name] = np.array(null_distributions[metric_name])
    
    # Compute p-values and effect sizes
    p_values = {}
    effect_sizes = {}
    
    for metric_name in metric_types:
        real_value = real_metrics.get(metric_name, np.nan)
        null_dist = null_distributions.get(metric_name, np.array([]))
        
        if len(null_dist) == 0 or np.isnan(real_value):
            p_values[metric_name] = np.nan
            effect_sizes[metric_name] = np.nan
            continue
        
        # P-value: fraction of permutations with value >= real value
        # (or <= for metrics where lower is better)
        if metric_name in ['davies_bouldin']:  # Lower is better
            p_value = np.mean(null_dist <= real_value)
        else:  # Higher is better (silhouette, etc.)
            p_value = np.mean(null_dist >= real_value)
        
        p_values[metric_name] = float(p_value)
        
        # Cohen's d effect size
        null_mean = float(null_dist.mean())
        null_std = float(null_dist.std())
        if null_std > 0:
            cohens_d = (real_value - null_mean) / null_std
        else:
            cohens_d = np.nan
        
        effect_sizes[metric_name] = float(cohens_d)
    
    # Print results
    console.print("\n[bold]Permutation Test Results:[/bold]")
    for metric_name in metric_types:
        real_val = real_metrics.get(metric_name, np.nan)
        null_dist_array = null_distributions[metric_name]
        
        if len(null_dist_array) > 0:
            console.print(f"\n  {metric_name}:")
            console.print(f"    Real value: {real_val:.4f}")
            console.print(f"    Null mean: {null_dist_array.mean():.4f} ± {null_dist_array.std():.4f}")
            console.print(f"    p-value: {p_values[metric_name]:.4f}")
            console.print(f"    Effect size (Cohen's d): {effect_sizes[metric_name]:.4f}")
            
            # Interpret effect size
            d = abs(effect_sizes[metric_name])
            if d < 0.2:
                effect_str = "negligible"
            elif d < 0.5:
                effect_str = "small"
            elif d < 0.8:
                effect_str = "medium"
            else:
                effect_str = "large"
            console.print(f"    Effect interpretation: [green]{effect_str}[/green]")
    
    return {
        'real_metrics': real_metrics,
        'null_distributions': null_distributions,
        'p_values': p_values,
        'effect_sizes': effect_sizes,
        'n_permutations': n_permutations,
    }


def interpret_significance(p_value, alpha=0.05):
    """Interpret p-value for significance testing.
    
    Parameters:
        p_value: p-value from permutation test
        alpha: Significance level (default 0.05)
        
    Returns:
        str: Significance interpretation
    """
    if p_value < alpha / 1000:
        return "*** (p < 0.001)"
    elif p_value < alpha / 100:
        return "** (p < 0.01)"
    elif p_value < alpha:
        return "* (p < 0.05)"
    else:
        return "ns (not significant)"


__all__ = [
    'compute_clustering_metrics',
    'permute_labels',
    'run_permutation_tests',
    'interpret_significance',
]

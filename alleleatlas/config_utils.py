"""Configuration validation and helper utilities for the refactored pipeline."""

from pathlib import Path
from rich.console import Console

console = Console()


def validate_config(config):
    """Validate ClusteringConfig object.
    
    Parameters:
        config: ClusteringConfig instance
        
    Raises:
        ValueError: If configuration is invalid
    """
    errors = []
    
    # Check required fields
    if not config.cgmlst_profiles:
        errors.append("cgmlst_profiles is required")
    elif not Path(config.cgmlst_profiles).exists():
        errors.append(f"Input profile file not found: {config.cgmlst_profiles}")
    
    if not config.outdir:
        errors.append("outdir is required")
    
    # Check numeric parameters
    if config.dbscan_min_samples < 1:
        errors.append("dbscan_min_samples must be >= 1")
    
    if config.dbscan_eps is not None and config.dbscan_eps < 0:
        errors.append("dbscan_eps must be >= 0")
    
    if config.min_cluster_size < 1:
        errors.append("min_cluster_size must be >= 1")
    
    if config.n_permutations < 10:
        errors.append("n_permutations should be >= 10 for meaningful statistics")
    
    if config.nproc < 1:
        errors.append("nproc must be >= 1")
    
    # Check layout method
    valid_layouts = ['spring', 'spectral', 'equidistant', 'circular']
    if config.layout_method not in valid_layouts:
        errors.append(f"layout_method must be one of: {valid_layouts}")
    
    # Check UPGMA method
    valid_methods = ['single', 'complete', 'average', 'weighted', 'centroid']
    if config.upgma_method not in valid_methods:
        errors.append(f"upgma_method must be one of: {valid_methods}")
    
    if errors:
        console.print("[bold red]Configuration Errors:[/bold red]")
        for error in errors:
            console.print(f"  ✗ {error}")
        raise ValueError(f"Invalid configuration: {len(errors)} error(s)")
    
    console.print("[green]✓[/green] Configuration validated")


def print_config_summary(config):
    """Print a summary of the configuration."""
    console.print("\n[bold]Configuration Summary:[/bold]")
    console.print(f"  Input: {config.cgmlst_profiles}")
    console.print(f"  Output: {config.outdir}")
    console.print("\n  [bold]DBSCAN:[/bold]")
    console.print(f"    eps: {config.dbscan_eps if config.dbscan_eps else 'auto-detect'}")
    console.print(f"    min_samples: {config.dbscan_min_samples}")
    console.print("\n  [bold]UPGMA:[/bold]")
    console.print(f"    method: {config.upgma_method}")
    console.print("\n  [bold]Permutation Testing:[/bold]")
    console.print(f"    enabled: {config.run_permutation_tests_flag}")
    console.print(f"    n_permutations: {config.n_permutations}")
    console.print("\n  [bold]Network Visualization:[/bold]")
    console.print(f"    layout: {config.layout_method}")
    console.print(f"    edge_threshold: {config.edge_threshold if config.edge_threshold else 'median'}")
    console.print(f"    edge_percentile: {config.edge_percentile if config.edge_percentile else 'none'}")
    console.print("\n  [bold]General:[/bold]")
    console.print(f"    processes: {config.nproc}")
    console.print(f"    force_recompute: {config.force_recompute}")


__all__ = [
    'validate_config',
    'print_config_summary',
]

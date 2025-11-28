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
    if config.k < 1:
        errors.append("k must be >= 1")
    valid_backends = ['exact', 'usearch']
    if config.backend not in valid_backends:
        errors.append(f"backend must be one of: {valid_backends}")
    if config.nproc < 1:
        errors.append("nproc must be >= 1")
    
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
    console.print(f"  Backend: {config.backend}")
    console.print(f"  k neighbors: {config.k}")
    console.print("\n  [bold]General:[/bold]")
    console.print(f"    processes: {config.nproc}")
    console.print(f"    force_recompute: {config.force_recompute}")


__all__ = [
    'validate_config',
    'print_config_summary',
]

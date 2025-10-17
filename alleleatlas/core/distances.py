"""Distance matrix computation for cgMLST profiles.

Handles parallel computation of pairwise and dual distance matrices
from normalized allelic profiles.
"""

import multiprocessing as mp
from multiprocessing import Pool
from rich.console import Console

# Force 'fork' mode on macOS for better compatibility with SharedArray and numba
if mp.get_start_method() != 'fork':
    try:
        mp.set_start_method('fork', force=True)
    except RuntimeError:
        # Already set, ignore
        pass

from alleleatlas.cluster.pHierCC import prepare_mat
from alleleatlas.cluster.getDistance import getDistance

console = Console()

# Threshold for switching to batch mode (memory-efficient)
LARGE_DATASET_THRESHOLD = 10000  # profiles


def compute_distance_matrices(normalized_path, nproc=4, allowed_missing=0.05):
    """Compute both dual and pairwise distance matrices from normalized profile.
    
    Parameters:
        normalized_path (str): Path to normalized cgMLST profile TSV
        nproc (int): Number of parallel processes (default: 4, max 8 for large datasets)
        
    Returns:
        tuple: (dist_dual, dist_p) - Distance matrices for clustering and evaluation
        
    Notes:
        - dist_dual: Used for hierarchical clustering (includes count information)
        - dist_p: Used for evaluation metrics (pairwise distances only)
        - Uses multiprocessing for faster computation on large profiles
        - Allowed missing loci: 3% of total loci
        - For very large datasets (>10k profiles), uses memory-efficient batch mode
    """
    console.print('Computing distance matrix...')
    
    mat, names = prepare_mat(normalized_path)
    n_profiles = mat.shape[0]
    
    # Limit number of processes for large datasets to prevent semaphore exhaustion
    actual_nproc = min(nproc, 4 if n_profiles > LARGE_DATASET_THRESHOLD else 8)
    if actual_nproc != nproc:
        console.print(f'  [yellow]Large dataset detected ({n_profiles} profiles)')
        console.print(f'  Reducing processes from {nproc} to {actual_nproc} to manage memory[/yellow]')
    
    pool = Pool(actual_nproc)
    try:
        dist_dual = getDistance(mat, 'dual_dist', pool, start=0, allowed_missing=allowed_missing)
        dist_p = getDistance(mat, 'p_dist', pool, start=0, allowed_missing=allowed_missing)
    finally:
        pool.close()
        pool.join()
    
    console.print('  ✓ Distance matrices computed')
    return dist_dual, dist_p, mat, names


# Backward compatibility alias
_compute_distances = compute_distance_matrices


__all__ = ['compute_distance_matrices']

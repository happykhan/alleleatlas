"""Distance matrix computation for cgMLST profiles.

Handles parallel computation of pairwise and dual distance matrices
from normalized allelic profiles.
"""

from multiprocessing import Pool
from rich.console import Console

from alleleatlas.cluster.pHierCC import prepare_mat
from alleleatlas.cluster.getDistance import getDistance

console = Console()


def compute_distance_matrices(normalized_path, nproc=4):
    """Compute both dual and pairwise distance matrices from normalized profile.
    
    Parameters:
        normalized_path (str): Path to normalized cgMLST profile TSV
        nproc (int): Number of parallel processes (default: 4)
        
    Returns:
        tuple: (dist_dual, dist_p) - Distance matrices for clustering and evaluation
        
    Notes:
        - dist_dual: Used for hierarchical clustering (includes count information)
        - dist_p: Used for evaluation metrics (pairwise distances only)
        - Uses multiprocessing for faster computation on large profiles
        - Allowed missing loci: 3% of total loci
    """
    console.print('Computing distance matrix...')
    
    mat, names = prepare_mat(normalized_path)
    pool = Pool(nproc)
    try:
        dist_dual = getDistance(mat, 'dual_dist', pool, start=0, allowed_missing=0.03)
        dist_p = getDistance(mat, 'p_dist', pool, start=0, allowed_missing=0.03)
    finally:
        pool.close()
        pool.join()
    
    console.print('  ✓ Distance matrices computed')
    return dist_dual, dist_p


# Backward compatibility alias
_compute_distances = compute_distance_matrices


__all__ = ['compute_distance_matrices']

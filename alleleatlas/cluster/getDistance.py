"""
Optimized distance matrix computation with numba parallelization and shared memory.

Improvements over baseline:
1. Numba parallel loops (prange) for dual_dist and p_dist
2. Precomputed masks and counts to avoid repeated checks
3. Larger chunking for reduced IPC overhead
4. Compact dtypes (int16) for reduced memory bandwidth
5. Caching of numba compiled functions

Note: The current implementation uses starmap for simplicity. For further optimization,
consider using Pool(initializer=...) pattern to attach shared arrays once per worker.
"""
from tempfile import NamedTemporaryFile
import numpy as np
import numba as nb
from numba import prange
import SharedArray as sa


def getDistance(data, func_name, pool, start=0, allowed_missing=0.0):
    """
    Compute distance matrix using parallel numba functions.
    
    Parameters:
        data: ndarray of allele profiles (N_profiles x N_loci)
        func_name: 'dual_dist' or 'p_dist'
        pool: multiprocessing.Pool instance
        start: row index to start from
        allowed_missing: fraction of missing loci allowed
    
    Returns:
        Distance matrix array (N-start x N x 2)
    """
    with NamedTemporaryFile(dir='.', prefix='HCC_') as file:
        prefix = 'file://{0}'.format(file.name)
        
        # Create shared arrays
        mat_buf = '{0}.mat.sa'.format(prefix)
        mat = sa.create(mat_buf, shape=data.shape, dtype=np.int16)  # Use int16 to save bandwidth
        mat[:] = data[:]
        
        # Precompute masks and counts in main process (only once)
        present_mask_buf = '{0}.mask.sa'.format(prefix)
        present_mask = sa.create(present_mask_buf, shape=data.shape, dtype=np.uint8)
        present_mask[:] = (data > 0).astype(np.uint8)
        
        present_count_buf = '{0}.count.sa'.format(prefix)
        present_count = sa.create(present_count_buf, shape=(data.shape[0],), dtype=np.int32)
        present_count[:] = np.sum(data > 0, axis=1).astype(np.int32)
        
        dist_buf = '{0}.dist.sa'.format(prefix)
        dist = sa.create(dist_buf, shape=[mat.shape[0] - start, mat.shape[0], 2], dtype=np.int32)
        dist[:] = 0
        
        # Run parallel computation
        __parallel_dist(
            mat_buf, func_name, dist_buf, present_mask_buf, present_count_buf,
            mat.shape, pool, start, allowed_missing
        )
        
        # Copy result before cleanup
        result = np.array(dist)
        
        # Cleanup
        sa.delete(mat_buf)
        sa.delete(dist_buf)
        sa.delete(present_mask_buf)
        sa.delete(present_count_buf)
    
    return result


def _worker_task(mat_buf, mask_buf, count_buf, dist_buf, func_name, s, e, start, allowed_missing):
    """
    Worker task: attach shared arrays and compute distances for a range.
    
    This is called per task by pool.starmap. For ultimate optimization,
    use Pool(initializer=) to attach once per worker instead.
    """
    # Attach shared arrays in this worker
    mat = sa.attach(mat_buf)
    dist = sa.attach(dist_buf)
    masks = sa.attach(mask_buf)
    counts = sa.attach(count_buf)
    
    try:
        # Get the appropriate distance function
        if func_name == 'dual_dist':
            func = dual_dist_parallel
        elif func_name == 'p_dist':
            func = p_dist_parallel
        else:
            raise ValueError(f"Unknown distance function: {func_name}")
        
        # Compute distances (skip first column which is index/strain)
        mat_view = mat[:, 1:]
        mask_view = masks[:, 1:]
        d = func(mat_view, mask_view, counts, s, e, allowed_missing)
        
        # Write result to shared distance array
        dist[(s - start):(e - start)] = d
    finally:
        # Detach (not strictly necessary but good practice)
        del mat, dist, masks, counts


def __parallel_dist(mat_buf, func_name, dist_buf, mask_buf, count_buf,
                    mat_shape, pool, start=0, allowed_missing=0.0):
    """
    Distribute work across pool with reduced IPC overhead.
    
    Uses larger chunks (1 per worker instead of many tiny tasks) and
    avoids repeated attachment by calling task with all buffer names.
    """
    n_pool = len(pool._pool)
    
    # Create balanced chunks covering [start, n_profiles)
    n_profiles = mat_shape[0]
    chunks = _make_chunks(n_profiles, n_pool, start=start)
    
    # Create task list with all buffer references
    tasks = [
        (mat_buf, mask_buf, count_buf, dist_buf, func_name, s, e, start, allowed_missing)
        for s, e in chunks
    ]
    
    # Run tasks in parallel
    pool.starmap(_worker_task, tasks)
    return


def _make_chunks(n, n_workers, start=0):
    """
    Create balanced (s, e) ranges for parallel work.
    
    Parameters:
        n: total number of items
        n_workers: number of workers
        start: starting index
    
    Returns:
        List of (s, e) tuples covering [start, n)
    """
    total = n - start
    base = total // n_workers
    extras = total % n_workers
    ranges = []
    s = start
    for i in range(n_workers):
        e = s + base + (1 if i < extras else 0)
        if e <= s:
            e = s + 1
        ranges.append((s, e))
        s = e
    return ranges


@nb.njit(parallel=True, fastmath=True, cache=True)
def dual_dist_parallel(mat, present_mask, present_count, s, e, allowed_missing):
    """
    Parallel dual distance computation using numba prange.
    
    Optimizations:
    - Uses precomputed masks and counts (no repeated checks)
    - prange for outer loop parallelization
    - fastmath=True for speed (safe for this application)
    - cache=True to reuse compiled code
    
    Parameters:
        mat: int16 array (N_profiles x N_loci)
        present_mask: uint8 array (N_profiles x N_loci) 1=present, 0=missing
        present_count: int32 array (N_profiles,) count of present per row
        s, e: range of rows to compute
        allowed_missing: fraction of missing loci allowed
    
    Returns:
        Distance array (e-s x N_profiles x 2)
    """
    N = mat.shape[0]
    L = mat.shape[1]
    out = np.zeros((e - s, N, 2), dtype=np.int32)
    
    # Precompute floating-point constants to avoid repeated creation
    am = np.float64(allowed_missing)
    eps = np.float64(1e-12)
    L_float = np.float64(L)
    
    for ii in prange(s, e):
        i_idx = ii - s
        ql = np.float64(present_count[ii])
        
        # Iterate j from 0..ii-1
        for j in range(ii):
            rl = np.float64(present_count[j])
            ad = np.float64(0.0)
            al = eps
            
            # Count agreements/disagreements
            for k in range(L):
                if present_mask[j, k]:
                    # mat[j, k] is present
                    if present_mask[ii, k]:
                        # Both present
                        al += np.float64(1.0)
                        if mat[ii, k] != mat[j, k]:
                            ad += np.float64(1.0)
            
            # Adjustment 1
            ll2 = ql - am * L_float
            if ll2 > al:
                ad += (ll2 - al)
                al = ll2
            out[i_idx, j, 1] = np.int32(ad / al * L_float + np.float64(0.5))
            
            # Adjustment 2
            ll = np.maximum(ql, rl) - am * L_float
            if ll > al:
                ad += (ll - al)
                al = ll
            out[i_idx, j, 0] = np.int32(ad / al * L_float + np.float64(0.5))
    
    return out


@nb.njit(parallel=True, fastmath=True, cache=True)
def p_dist_parallel(mat, present_mask, present_count, s, e, allowed_missing):
    """
    Parallel pairwise distance computation using numba prange.
    
    Optimizations:
    - Uses precomputed masks (no repeated mat[j,k]>0 checks)
    - prange for outer loop parallelization
    - fastmath=True for log computation speed
    - cache=True for compiled code reuse
    
    Parameters:
        mat: int16 array (N_profiles x N_loci)
        present_mask: uint8 array (N_profiles x N_loci)
        present_count: int32 array (N_profiles,) [unused but kept for API consistency]
        s, e: range of rows to compute
        allowed_missing: [unused but kept for API consistency]
    
    Returns:
        Distance array (e-s x N_profiles x 2)
    """
    N = mat.shape[0]
    L = mat.shape[1]
    out = np.zeros((e - s, N, 2), dtype=np.int32)
    
    L_float = np.float64(L)
    
    for ii in prange(s, e):
        i_idx = ii - s
        
        for j in range(ii):
            ad = np.float64(0.0)
            al = np.float64(0.0)
            
            # Count agreements/disagreements using precomputed mask
            for k in range(L):
                if present_mask[j, k] and present_mask[ii, k]:
                    al += np.float64(1.0)
                    if mat[ii, k] != mat[j, k]:
                        ad += np.float64(1.0)
            
            # Compute log-likelihood distance
            denom = al + np.float64(1.0)
            tmp = (ad + np.float64(0.5)) / denom
            
            # Guard against log(0)
            if tmp >= np.float64(1.0):
                tmp = np.float64(1.0) - np.float64(1e-12)
            
            val = -np.log(np.float64(1.0) - tmp) * L_float * np.float64(100.0) + np.float64(0.5)
            out[i_idx, j, 0] = np.int32(val)
    
    return out


# Legacy non-parallel versions for testing/fallback
@nb.jit(nopython=True)
def dual_dist(mat, s, e, allowed_missing=0.05):
    """Legacy serial implementation for reference/fallback."""
    dist = np.zeros((e - s, mat.shape[0], 2), dtype=np.int32)
    n_loci = mat.shape[1]
    for i in range(s, e):
        ql = np.sum(mat[i] > 0)
        for j in range(i):
            rl, ad, al = 0., 1e-4, 1e-4
            for k in range(n_loci):
                if mat[j, k] > 0:
                    rl += 1
                    if mat[i, k] > 0:
                        al += 1
                        if mat[i, k] != mat[j, k]:
                            ad += 1
            ll = max(ql, rl) - allowed_missing * n_loci
            ll2 = ql - allowed_missing * n_loci

            if ll2 > al:
                ad += ll2 - al
                al = ll2
            dist[i - s, j, 1] = int(ad / al * n_loci + 0.5)

            if ll > al:
                ad += ll - al
                al = ll
            dist[i - s, j, 0] = int(ad / al * n_loci + 0.5)
    return dist


@nb.jit(nopython=True)
def p_dist(mat, s, e, allowed_missing=0.05):
    """Legacy serial implementation for reference/fallback."""
    dist = np.zeros((e - s, mat.shape[0], 2), dtype=np.int32)
    n_loci = mat.shape[1]
    for i in range(s, e):
        for j in range(i):
            ad, al = 0., 0.
            for k in range(n_loci):
                if mat[j, k] > 0:
                    if mat[i, k] > 0:
                        al += 1
                        if mat[i, k] != mat[j, k]:
                            ad += 1
            dist[i - s, j, 0] = int(-np.log(1. - (ad + 0.5) / (al + 1.0)) * n_loci * 100. + 0.5)
    return dist

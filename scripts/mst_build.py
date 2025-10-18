import numpy as np
import networkx as nx
from typing import Optional
from rich import print as rprint
import time
import os
from rich.console import Console
console = Console()

def get_memory_usage():
    """Get current process memory usage in MB (platform-specific)."""
    try:
        import psutil
        process = psutil.Process(os.getpid())
        return process.memory_info().rss / 1024 / 1024
    except ImportError:
        # Fallback: estimate via resource module
        try:
            import resource
            return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        except Exception:
            return 0  # Can't track memory on this platform


def generate_synthetic_distance_matrix(n_nodes: int, seed: int = 42) -> np.ndarray:
    """Generate a symmetric distance matrix for benchmarking.
    
    Creates a random symmetric matrix with values between 0 and 100,
    representing pairwise distances between samples.
    """
    np.random.seed(seed)
    
    # Generate random upper triangle
    upper = np.random.randint(0, 101, size=(n_nodes, n_nodes))
    
    # Make symmetric
    dist_matrix = np.maximum(upper, upper.T).astype(np.float32)
    
    # Zero diagonal
    np.fill_diagonal(dist_matrix, 0)
    
    return dist_matrix

def build_mst_from_distances_alt(
    distance_matrix: np.ndarray,
    *,
    use_knn: Optional[int] = None,
    knn_algorithm: str = "auto",
    dtype: np.dtype = np.float32,
    verbose: bool = True,
) -> nx.Graph:
    """
    Build MST from a symmetric distance matrix.

    Parameters
    ----------
    distance_matrix :
        2D square numpy array (n, n). Must be symmetric, diagonal can be 0 or np.nan.
    use_knn :
        If set to an int k, build a k-nearest-neighbor sparse graph first (k neighbors per node)
        and compute MST on that sparse graph. This reduces memory and time dramatically for large n.
        If None, use the full pairwise edges (but stored as sparse COO to avoid a huge NetworkX graph).
    knn_algorithm :
        Passed to sklearn.neighbors.NearestNeighbors (if available): 'auto', 'ball_tree', 'kd_tree', 'brute'
    dtype :
        Use float32 if precision is acceptable.
    verbose :
        Whether to print progress.

    Returns
    -------
    G : networkx.Graph
        MST with attributes 'weight' and 'distance' on edges.
    """
    import scipy.sparse as sp
    from scipy.sparse.csgraph import minimum_spanning_tree

    n = distance_matrix.shape[0]
    if verbose:
        rprint(f"[cyan]Building MST from {n}x{n} distance matrix...[/cyan]")

    # Cast to desired dtype and ensure C contiguous
    if distance_matrix.dtype != dtype:
        distance_matrix = distance_matrix.astype(dtype, copy=False)
    else:
        distance_matrix = np.ascontiguousarray(distance_matrix)

    # Option 1: k-NN pruning (recommended for large n)
    if use_knn is not None and use_knn > 0 and use_knn < n:
        try:
            from sklearn.neighbors import NearestNeighbors
            if verbose:
                rprint(f"  Using exact k-NN graph with k={use_knn} (sklearn)...")
            # Use n_jobs=-1 to parallelize neighbor search
            nn = NearestNeighbors(n_neighbors=min(use_knn + 1, n), algorithm=knn_algorithm, n_jobs=-1)
            # For NearestNeighbors we provide the matrix of features; but we have distances already.
            # We'll instead use a trick: for each row find the k smallest distances (excluding self).
            # This is done in pure numpy below (vectorized partial sort per row):
            # We assume matrix is dense and fits memory.
            # Mask self distances to +inf to exclude them
            mat = distance_matrix.copy()
            np.fill_diagonal(mat, np.inf)
            # Use argpartition to find k smallest per row (fast)
            kplus = min(use_knn, n - 1)
            idx_k = np.argpartition(mat, kth=kplus, axis=1)[:, :kplus]  # shape (n, k)
            rows = np.repeat(np.arange(n), kplus)
            cols = idx_k.reshape(-1)
            data = mat[np.arange(n)[:, None], idx_k].reshape(-1)
            # Filter out invalid edges (NaN or negative)
            valid_mask = (~np.isnan(data)) & (data >= 0) & (data != np.inf)
            rows = rows[valid_mask]
            cols = cols[valid_mask]
            data = data[valid_mask]

            # Make symmetric: if A->B exists, ensure B->A exists too by adding the reversed edges
            rows_sym = np.concatenate([rows, cols])
            cols_sym = np.concatenate([cols, rows])
            data_sym = np.concatenate([data, data])
            coo = sp.coo_matrix((data_sym, (rows_sym, cols_sym)), shape=(n, n))
            csr = coo.tocsr()
            if verbose:
                rprint(f"    k-NN edges: {coo.nnz // 2:,} unique undirected edges (stored doubled).")
        except Exception as e:
            # sklearn not available or failed; fallback to full-sparse approach
            if verbose:
                rprint("[yellow]sklearn NearestNeighbors unavailable or failed, falling back to full-sparse.[/yellow]")
            use_knn = None  # fall through to full approach

    # Option 2: full sparse representation (vectorized)
    if use_knn is None or use_knn <= 0:
        if verbose:
            rprint("  Building sparse COO from upper triangle (vectorized)...")
        # vectorized extraction of upper-triangle indices (i<j)
        iu = np.triu_indices(n, k=1)
        data_full = distance_matrix[iu]
        # mask valid edges (non-nan, non-negative, maybe non-inf)
        valid = (~np.isnan(data_full)) & (data_full >= 0) & (np.isfinite(data_full))
        rows = iu[0][valid]
        cols = iu[1][valid]
        data = data_full[valid]
        if verbose:
            rprint(f"    Valid edges found: {data.size:,}")

        # Build symmetric COO (both directions), required by csgraph routines
        rows_sym = np.concatenate([rows, cols])
        cols_sym = np.concatenate([cols, rows])
        data_sym = np.concatenate([data, data])
        coo = sp.coo_matrix((data_sym, (rows_sym, cols_sym)), shape=(n, n))
        csr = coo.tocsr()

    # Compute minimum spanning tree on sparse CSR (Prim/Kruskal internal)
    if verbose:
        rprint("  Computing minimum spanning tree (scipy.sparse.csgraph.minimum_spanning_tree)...")
    mst_sparse = minimum_spanning_tree(csr)
    # mst_sparse is sparse in CSR with directed edges (i->j). Convert to COO for iteration.
    mst_coo = mst_sparse.tocoo()

    # Build small NetworkX graph from MST edges (only ~n-1 edges)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    # The MST could have directed entries; we'll treat them as undirected unique edges
    # Iterate and add undirected edge only once
    seen = set()
    for i, j, w in zip(mst_coo.row, mst_coo.col, mst_coo.data):
        a, b = int(i), int(j)
        if a == b:
            continue
        key = (min(a, b), max(a, b))
        if key in seen:
            continue
        seen.add(key)
        G.add_edge(key[0], key[1], weight=float(w), distance=float(w))
    # Add node attributes
    for node in G.nodes():
        G.nodes[node]["name"] = f"S{node}"
        G.nodes[node]["count"] = 1
    # Mark MST edges
    for u, v in G.edges():
        G[u][v]["is_mst"] = True

    if verbose:
        rprint(f"[green]Built MST[/green]: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G




# -----------------------
# MST construction
# -----------------------
def build_mst_from_distances(distance_matrix: np.ndarray) -> nx.Graph:
    """Create MST from distance matrix using optimized approach for large datasets.
    
    For large datasets (>10k profiles), uses scipy.sparse_coo for efficient
    edge representation and Prim's algorithm via NetworkX.
    """
    n = distance_matrix.shape[0]
    console.print(f"  Building MST from {n}x{n} distance matrix...")
    
    # For large datasets, use sparse representation to avoid creating complete graph
    if n > 5000:
        console.print("    Using sparse edge representation for large dataset...")
        from scipy.sparse import coo_matrix
        
        # Convert distance matrix to COO format, keeping only valid edges
        # This avoids storing all n² edges in memory
        rows = []
        cols = []
        data = []
        
        for i in range(n):
            # Only check upper triangle to save iterations
            for j in range(i + 1, n):
                d = distance_matrix[i, j]
                # Skip NaN and negative distances
                if d >= 0 and not np.isnan(d):
                    rows.append(i)
                    cols.append(j)
                    data.append(d)
        
        console.print(f"    Edges to process: {len(data):,}")
        
        # Create sparse COO matrix
        coo = coo_matrix((data, (rows, cols)), shape=(n, n))
        
        # Convert to CSR for efficient MST computation
        csr = coo.tocsr()
        
        # Use scipy's minimum spanning tree for sparse graphs
        from scipy.sparse.csgraph import minimum_spanning_tree
        mst_sparse = minimum_spanning_tree(csr)
        
        # Convert back to NetworkX graph
        G = nx.Graph()
        G.add_nodes_from(range(n))
        
        # Extract edges from sparse MST
        mst_coo = mst_sparse.tocoo()
        for i, j, w in zip(mst_coo.row, mst_coo.col, mst_coo.data):
            if i < j:  # Only add each edge once (graph is undirected)
                G.add_edge(i, j, weight=float(w), distance=float(w))
    
    else:
        # For smaller datasets, use the simpler NetworkX approach
        G = nx.Graph()
        G.add_nodes_from(range(n))
        
        console.print("  Adding edges to the graph...")
        edge_count = 0
        for i in range(n):
            for j in range(i + 1, n):
                d = distance_matrix[i, j]
                if d >= 0 and not np.isnan(d):
                    G.add_edge(i, j, weight=float(d), distance=float(d))
                    edge_count += 1
        console.print(f"    Added {edge_count:,} edges")
        
        console.print("[bold]Computing MST with Prim's algorithm...[/bold]")
        G = nx.minimum_spanning_tree(G, weight="weight")
    
    # Add node attributes
    for node in G.nodes():
        G.nodes[node]["name"] = f"S{node}"
        G.nodes[node]["count"] = 1
    
    # Mark MST edges (all edges in result are MST edges)
    for u, v in G.edges():
        G[u][v]["is_mst"] = True
    
    rprint(f"[bold]Built MST[/bold]: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


# -----------------------
# Benchmarking
# -----------------------
def benchmark_method(method_name: str, method_func, distance_matrix: np.ndarray) -> dict:
    """Benchmark a single MST construction method."""
    n = distance_matrix.shape[0]
    
    mem_before = get_memory_usage()
    time_start = time.perf_counter()
    
    try:
        G = method_func(distance_matrix)
        time_elapsed = time.perf_counter() - time_start
        mem_after = get_memory_usage()
        mem_used = mem_after - mem_before
        
        return {
            'method': method_name,
            'nodes': n,
            'time_sec': time_elapsed,
            'memory_mb': mem_used,
            'edges': G.number_of_edges(),
            'status': 'success'
        }
    except Exception as e:
        time_elapsed = time.perf_counter() - time_start
        return {
            'method': method_name,
            'nodes': n,
            'time_sec': time_elapsed,
            'memory_mb': 0,
            'edges': 0,
            'status': f'failed: {str(e)[:50]}'
        }


def method_1_full_sparse(dist_matrix: np.ndarray) -> nx.Graph:
    """Method 1: Full sparse approach (vectorized)."""
    n = dist_matrix.shape[0]
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import minimum_spanning_tree
    
    # Vectorized extraction of upper triangle
    iu = np.triu_indices(n, k=1)
    data_full = dist_matrix[iu]
    valid = (~np.isnan(data_full)) & (data_full >= 0) & (np.isfinite(data_full))
    rows = iu[0][valid]
    cols = iu[1][valid]
    data = data_full[valid]
    
    # Make symmetric
    rows_sym = np.concatenate([rows, cols])
    cols_sym = np.concatenate([cols, rows])
    data_sym = np.concatenate([data, data])
    coo = coo_matrix((data_sym, (rows_sym, cols_sym)), shape=(n, n))
    csr = coo.tocsr()
    
    # Compute MST
    mst_sparse = minimum_spanning_tree(csr)
    mst_coo = mst_sparse.tocoo()
    
    # Build NetworkX graph
    G = nx.Graph()
    G.add_nodes_from(range(n))
    seen = set()
    for i, j, w in zip(mst_coo.row, mst_coo.col, mst_coo.data):
        a, b = int(i), int(j)
        if a == b:
            continue
        key = (min(a, b), max(a, b))
        if key in seen:
            continue
        seen.add(key)
        G.add_edge(key[0], key[1], weight=float(w), distance=float(w))
    
    return G


def method_2_knn_sparse(dist_matrix: np.ndarray, k: int = 50) -> nx.Graph:
    """Method 2: k-NN pruning approach."""
    n = dist_matrix.shape[0]
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import minimum_spanning_tree
    
    # k-NN using argpartition
    mat = dist_matrix.copy()
    np.fill_diagonal(mat, np.inf)
    kplus = min(k, n - 1)
    idx_k = np.argpartition(mat, kth=kplus, axis=1)[:, :kplus]
    rows = np.repeat(np.arange(n), kplus)
    cols = idx_k.reshape(-1)
    data = mat[np.arange(n)[:, None], idx_k].reshape(-1)
    
    # Filter invalid
    valid_mask = (~np.isnan(data)) & (data >= 0) & (data != np.inf)
    rows = rows[valid_mask]
    cols = cols[valid_mask]
    data = data[valid_mask]
    
    # Make symmetric
    rows_sym = np.concatenate([rows, cols])
    cols_sym = np.concatenate([cols, rows])
    data_sym = np.concatenate([data, data])
    coo = coo_matrix((data_sym, (rows_sym, cols_sym)), shape=(n, n))
    csr = coo.tocsr()
    
    # Compute MST
    mst_sparse = minimum_spanning_tree(csr)
    mst_coo = mst_sparse.tocoo()
    
    # Build NetworkX graph
    G = nx.Graph()
    G.add_nodes_from(range(n))
    seen = set()
    for i, j, w in zip(mst_coo.row, mst_coo.col, mst_coo.data):
        a, b = int(i), int(j)
        if a == b:
            continue
        key = (min(a, b), max(a, b))
        if key in seen:
            continue
        seen.add(key)
        G.add_edge(key[0], key[1], weight=float(w), distance=float(w))
    
    return G


def method_3_networkx_complete(dist_matrix: np.ndarray) -> nx.Graph:
    """Method 3: Complete graph in NetworkX (baseline, slower for large n)."""
    n = dist_matrix.shape[0]
    
    G = nx.Graph()
    G.add_nodes_from(range(n))
    
    edge_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = dist_matrix[i, j]
            if d >= 0 and not np.isnan(d):
                G.add_edge(i, j, weight=float(d), distance=float(d))
                edge_count += 1
    
    # Compute MST
    G = nx.minimum_spanning_tree(G, weight="weight")
    return G


def run_benchmarks():
    """Run comprehensive benchmarks across different dataset sizes."""
    import time
    
    rprint("\n[bold cyan]MST Construction Method Benchmarks[/bold cyan]")
    rprint("=" * 80)
    
    # Dataset sizes to test
    dataset_sizes = [100, 500, 1000, 2000, 5000, 10000, 20000]
    
    results = []
    
    for n_nodes in dataset_sizes:
        rprint(f"\n[bold]Testing with {n_nodes} nodes...[/bold]")
        dist_matrix = generate_synthetic_distance_matrix(n_nodes)
        
        # Test Method 1: Full Sparse
        rprint("  Method 1 (full sparse)...", end="")
        result1 = benchmark_method("Full Sparse", method_1_full_sparse, dist_matrix)
        rprint(f" {result1['time_sec']:.3f}s ({result1['memory_mb']:.1f}MB)")
        results.append(result1)
        
        # Test Method 2: k-NN (k=50)
        rprint("  Method 2 (k-NN, k=50)...", end="")
        result2 = benchmark_method("k-NN Sparse", lambda m: method_2_knn_sparse(m, k=50), dist_matrix)
        rprint(f" {result2['time_sec']:.3f}s ({result2['memory_mb']:.1f}MB)")
        results.append(result2)
        
        # Test Method 3: NetworkX (only for smaller sizes due to O(n²) complexity)
        if n_nodes <= 2000:
            rprint("  Method 3 (NetworkX complete)...", end="")
            result3 = benchmark_method("NetworkX Complete", method_3_networkx_complete, dist_matrix)
            rprint(f" {result3['time_sec']:.3f}s ({result3['memory_mb']:.1f}MB)")
            results.append(result3)
    
    # Print summary table
    rprint("\n[bold]Summary[/bold]")
    rprint("=" * 80)
    rprint(f"{'Method':<20} {'Nodes':<10} {'Time (s)':<12} {'Memory (MB)':<15} {'Status':<20}")
    rprint("-" * 80)
    
    for r in results:
        rprint(f"{r['method']:<20} {r['nodes']:<10} {r['time_sec']:<12.4f} {r['memory_mb']:<15.1f} {r['status']:<20}")
    
    # Analysis
    rprint("\n[bold]Analysis[/bold]")
    rprint("=" * 80)
    
    # Group by size
    from collections import defaultdict
    by_size = defaultdict(list)
    for r in results:
        by_size[r['nodes']].append(r)
    
    for size in sorted(by_size.keys()):
        group = by_size[size]
        successful = [r for r in group if r['status'] == 'success']
        if len(successful) > 1:
            fastest = min(successful, key=lambda r: r['time_sec'])
            slowest = max(successful, key=lambda r: r['time_sec'])
            speedup = slowest['time_sec'] / fastest['time_sec']
            rprint(f"n={size}: Fastest={fastest['method']} ({fastest['time_sec']:.3f}s), "
                   f"Slowest={slowest['method']} ({slowest['time_sec']:.3f}s), "
                   f"Speedup={speedup:.2f}x")


if __name__ == "__main__":
    run_benchmarks()

"""MST (Minimum Spanning Tree) visualization with sample counts and group coloring.

Generates network visualizations showing:
- Nodes: clusters/profiles with size proportional to sample count
- Edges: MST connections
- Colors: group/cluster assignments
"""

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx
from rich.console import Console

console = Console()


def build_mst_from_distances(distances, threshold=None):
    """Build Minimum Spanning Tree from distance matrix.
    
    Parameters:
        distances: scipy.sparse distance matrix or numpy array
        threshold: edge weight threshold (exclude edges above this)
        
    Returns:
        networkx.Graph: MST graph
    """
    # Convert to dense if sparse
    if hasattr(distances, 'toarray'):
        dist_array = distances.toarray()
    else:
        dist_array = np.asarray(distances)
    
    # Handle 3D distance format (dual_dist returns [dist, count] pairs)
    if len(dist_array.shape) == 3 and dist_array.shape[2] == 2:
        # dual_dist format: [..., 0] = distance, [..., 1] = count
        dist_array = dist_array[:, :, 0]
    
    # Build complete graph from distances
    n = dist_array.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(n))
    
    for i in range(n):
        for j in range(i + 1, n):
            weight = float(dist_array[i, j])  # Ensure scalar float
            if threshold is None or weight <= threshold:
                G.add_edge(i, j, weight=weight)
    
    # Get MST
    mst = nx.minimum_spanning_tree(G)
    return mst


def get_cluster_groups_from_hierarchy(hiercc_gz_path, threshold):
    """Extract cluster group assignments from HierCC clustering file.
    
    Parameters:
        hiercc_gz_path: path to hierCC.HierCC.gz file
        threshold: HC threshold (e.g., 400 for HC400)
        
    Returns:
        dict: mapping of sample_id → group_id
    """
    try:
        # Try to read from gzipped HierCC file
        hiercc_df = pd.read_csv(hiercc_gz_path, sep='\t', compression='gzip', low_memory=False)
        
        # Look for column matching threshold (e.g., 'HC400')
        threshold_col = f'HC{threshold}' if isinstance(threshold, int) else str(threshold)
        
        if threshold_col in hiercc_df.columns:
            groups = {}
            for idx, row in hiercc_df.iterrows():
                groups[idx] = row[threshold_col]
            return groups
    except Exception as e:
        console.print(f'[yellow]⚠[/yellow] Could not read groups from HierCC: {e}')
    
    # Fallback: no groups
    return {}


def draw_mst_network(
    mst,
    outpath,
    title,
    sample_counts=None,
    group_assignments=None,
    collapse_metadata=None,
    figsize=(12, 10)
):
    """Draw MST network with node sizes and group coloring.
    
    Parameters:
        mst: networkx.Graph (the MST)
        outpath: output PNG file path
        title: plot title
        sample_counts: dict mapping node_id → sample count
        group_assignments: dict mapping node_id → group/color
        collapse_metadata: metadata dict with sample_counts key
        figsize: figure size tuple
    """
    # Extract sample counts from metadata if available
    if sample_counts is None and collapse_metadata is not None:
        sample_counts = collapse_metadata.get('sample_counts', {})
    
    if sample_counts is None:
        sample_counts = {}
    
    # Node sizes scaled by sample count (or default 100)
    node_sizes = [
        max(100, sample_counts.get(node, 1) * 50)
        for node in mst.nodes()
    ]
    
    # Colors for groups - convert RGBA tuples to hex strings
    if group_assignments:
        # Create color mapping
        unique_groups = sorted(set(group_assignments.values()))
        num_colors = min(len(unique_groups), 20)
        colormap = plt.cm.get_cmap('tab20', num_colors)
        
        group_to_color_idx = {group: i % num_colors for i, group in enumerate(unique_groups)}
        node_color_indices = [
            group_to_color_idx.get(group_assignments.get(node, 0), 0)
            for node in mst.nodes()
        ]
        # Convert RGBA tuples to hex colors
        node_colors = [
            mcolors.to_hex(colormap(idx)) for idx in node_color_indices
        ]
    else:
        # Default blue color
        node_colors = 'steelblue'
    
    # Draw
    plt.figure(figsize=figsize)
    pos = nx.spring_layout(mst, k=2, iterations=50, seed=42)
    
    # Draw edges
    nx.draw_networkx_edges(mst, pos, alpha=0.3, width=1.0)
    
    # Draw nodes
    nx.draw_networkx_nodes(
        mst, pos,
        node_size=node_sizes,
        node_color=node_colors,
        alpha=0.9,
        edgecolors='black',
        linewidths=1.0
    )
    
    # Draw labels (optional - skip if too many nodes)
    if len(mst.nodes()) < 50:
        nx.draw_networkx_labels(mst, pos, font_size=8)
    
    plt.title(title, fontsize=14, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()


def run_mst_visualization(
    normalized_path,
    hiercc_path,
    outdir,
    collapse_threshold=0,
    line_threshold=None,
    collapse_metadata=None,
    precomputed_dist=None,
):
    """Generate MST visualization from clustering results.
    
    Parameters:
        normalized_path: path to normalized profile TSV
        hiercc_path: path to hierCC.HierCC.gz clustering file
        outdir: output directory
        collapse_threshold: HC threshold for node collapse (e.g., 400)
        line_threshold: HC threshold for MST edges (e.g., 1000)
        collapse_metadata: metadata from profile collapse operation
        precomputed_dist: precomputed distance matrix
        
    Returns:
        Path: path to output PNG file
    """
    outdir = Path(outdir)
    outpath = outdir / 'mst_visualization.png'
    
    # Load profile
    profile_df = pd.read_csv(normalized_path, sep='\t', low_memory=False)
    n_samples = len(profile_df)
    
    # Get sample counts
    sample_counts = {}
    if collapse_metadata and 'sample_counts' in collapse_metadata:
        sample_counts = collapse_metadata['sample_counts']
    else:
        # Default: each sample is 1
        sample_counts = {i: 1 for i in range(n_samples)}
    
    # Get group assignments from clustering
    group_assignments = {}
    if collapse_threshold is not None:
        try:
            group_assignments = get_cluster_groups_from_hierarchy(
                hiercc_path, collapse_threshold
            )
        except Exception as e:
            console.print(f'[yellow]⚠[/yellow] Could not load group assignments: {e}')
    
    # Build or use precomputed MST
    if precomputed_dist is not None:
        try:
            mst = build_mst_from_distances(precomputed_dist, threshold=line_threshold)
        except Exception as e:
            console.print(f'[yellow]⚠[/yellow] Error building MST from distances: {e}')
            return outpath
    else:
        console.print('[yellow]⚠[/yellow] No precomputed distances; MST visualization skipped')
        return outpath
    
    # Generate title
    title = f'MST Cluster Network (HC{collapse_threshold})'
    if line_threshold is not None:
        title += f' - Lines ≥{line_threshold}'
    title += f'\n{len(mst.nodes())} clusters, {len(mst.edges())} edges'
    
    # Draw
    draw_mst_network(
        mst,
        outpath,
        title,
        sample_counts=sample_counts,
        group_assignments=group_assignments,
        collapse_metadata=collapse_metadata
    )
    
    console.print(f'[green]✓[/green] MST visualization: {outpath}')
    return outpath


__all__ = ['run_mst_visualization', 'build_mst_from_distances', 'draw_mst_network']

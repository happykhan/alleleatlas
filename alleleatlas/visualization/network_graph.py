"""Network graph visualization with multiple layout options.

Provides NetworkX-based network visualization with flexible layout engines
(spring, spectral, equidistant) for exploring clustering results.
"""

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path
from rich.console import Console

console = Console()


class NetworkGraph:
    """Base class for network visualization."""
    
    def __init__(
        self,
        distance_matrix,
        labels=None,
        sample_names=None,
    ):
        """Initialize network graph.
        
        Parameters:
            distance_matrix: Symmetric distance matrix (n_samples x n_samples) or 3D array from getDistance
            labels: Cluster labels for coloring nodes
            sample_names: Optional sample names/IDs
            collapse_counts: Optional dict {sample_id: count} for node sizing
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
        
        self.distance_matrix = distance_matrix
        self.labels = labels
        self.sample_names = sample_names if sample_names is not None else [f"S{i}" for i in range(len(distance_matrix))]
        self.graph = None
        self.pos = None
    
    def build_graph(self, edge_threshold=None, edge_percentile=None):
        """Build NetworkX graph from distance matrix using Minimum Spanning Tree.
        
        Parameters:
            edge_threshold: Unused (MST always computed)
            edge_percentile: Unused (MST always computed)
        """
        console.print("Building network graph...")        
        # Create complete graph with distances as weights
        complete_graph = nx.Graph(self.distance_matrix)        
        # Get MST - this becomes our graph
        self.graph = nx.minimum_spanning_tree(complete_graph, weight='weight')
       
        return self.graph
    
    def compute_layout(self, method='tree', **kwargs):
        """Compute node positions using specified layout.
        
        Parameters:
            method: 'spring', 'spectral', 'circular', 'tree', 'radial', or 'hierarchical'
            use_mst_aware: Use MST-aware spring layout for better tree visualization
            **kwargs: Additional arguments for layout algorithm
        """
        if self.graph is None:
            raise ValueError("Build graph first with build_graph()")
        
        console.print(f"Computing {method} layout...")
        
        if method == 'spring':
            # Enhanced spring layout for MST - use edge weights inversely to respect tree structure
            self.pos = nx.spring_layout(
                self.graph,
                seed=42,
                iterations=100,
                weight='weight',  # Use edge weights
                **kwargs
            )

        
        elif method == 'tree':
            # Hierarchical tree layout (top-down)
            try:
                if nx.is_tree(self.graph):
                    # Find root (node with highest degree or center)
                    root = self._find_tree_root()
                    self.pos = self._tree_layout(root, level=0, width=1.0)
                else:
                    console.print("[yellow]Graph is not a tree, using spring layout[/yellow]")
                    self.pos = nx.spring_layout(self.graph, seed=42, iterations=100, k=2.0)
            except Exception as e:
                console.print(f"[yellow]Tree layout failed: {e}, using spring layout[/yellow]")
                self.pos = nx.spring_layout(self.graph, seed=42, iterations=100, k=2.0)
        
        elif method == 'radial':
            # Radial/circular tree layout (circular from center)
            try:
                if nx.is_tree(self.graph):
                    root = self._find_tree_root()
                    self.pos = self._radial_tree_layout(root)
                else:
                    console.print("[yellow]Graph is not a tree, using circular layout[/yellow]")
                    self.pos = nx.circular_layout(self.graph, **kwargs)
            except Exception as e:
                console.print(f"[yellow]Radial layout failed: {e}, using circular layout[/yellow]")
                self.pos = nx.circular_layout(self.graph, **kwargs)
        
        elif method == 'hierarchical':
            # Layered hierarchical layout (Sugiyama-style)
            try:
                self.pos = self._hierarchical_layout()
            except Exception as e:
                console.print(f"[yellow]Hierarchical layout failed: {e}, using spring layout[/yellow]")
                self.pos = nx.spring_layout(self.graph, seed=42, iterations=100, k=2.0)
        
        elif method == 'spectral':
            try:
                self.pos = nx.spectral_layout(self.graph, **kwargs)
            except Exception:
                console.print("[yellow]Spectral layout failed, falling back to spring[/yellow]")
                self.pos = nx.spring_layout(self.graph, seed=42, **kwargs)
        
        elif method == 'circular':
            self.pos = nx.circular_layout(self.graph, **kwargs)
        
        else:
            raise ValueError(f"Unknown layout method: {method}")
        
        console.print("  [green]✓[/green] Layout computed")
        return self.pos
    
    def _find_tree_root(self):
        """Find a good root node for tree layout (node with highest degree)."""
        if not self.graph or len(list(self.graph.nodes())) == 0:
            return 0
        return max(self.graph.nodes(), key=lambda x: self.graph.degree(x))  # type: ignore
    
    def _tree_layout(self, node, parent=None, level=0, pos=None, width=1.0, vert_gap=1.0, vert_loc=0.0, xcenter=0.0):
        """Compute positions for a tree layout (top-down)."""
        if pos is None:
            pos = {}
        
        pos[node] = (xcenter, vert_loc)
        neighbors = list(self.graph.neighbors(node))  # type: ignore
        
        if parent is not None:
            neighbors.remove(parent)
        
        if len(neighbors) > 0:
            dx = width / (len(neighbors) + 1)
            nextx = xcenter - width / 2 - dx / 2
            for neighbor in neighbors:
                nextx += dx
                self._tree_layout(neighbor, node, level + 1, pos, width=dx, 
                                 vert_gap=vert_gap, vert_loc=vert_loc - vert_gap, xcenter=nextx)
        
        return pos
    
    def _radial_tree_layout(self, root, level=0.0, angle=0.0, angle_width=2 * np.pi, distance=1.0, pos=None, parent=None):
        """Compute positions for radial/circular tree layout."""
        if pos is None:
            pos = {}
        
        # Place node at current angle and distance
        x = distance * np.cos(angle)
        y = distance * np.sin(angle)
        pos[root] = (x, y)
        
        neighbors = list(self.graph.neighbors(root))  # type: ignore
        if parent is not None and parent in neighbors:
            neighbors.remove(parent)
        
        if len(neighbors) > 0:
            # Divide angle space among children
            angle_per_child = angle_width / len(neighbors)
            start_angle = angle - angle_width / 2
            
            for i, child in enumerate(neighbors):
                child_angle = start_angle + (i + 0.5) * angle_per_child
                child_distance = distance + 1.0  # Increase distance at each level
                self._radial_tree_layout(child, level + 1, child_angle, angle_per_child,
                                        child_distance, pos, root)
        
        return pos
    
    def _hierarchical_layout(self):
        """Compute layered hierarchical layout."""
        pos = {}
        
        # Find a good root
        if self.graph and nx.is_tree(self.graph):  # type: ignore
            root = self._find_tree_root()
        else:
            root = 0
        
        # Assign levels via BFS
        visited = set()
        queue = [(root, 0)]
        levels = {}
        
        while queue:
            node, level = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            
            if level not in levels:
                levels[level] = []
            levels[level].append(node)
            
            for neighbor in self.graph.neighbors(node):  # type: ignore
                if neighbor not in visited:
                    queue.append((neighbor, level + 1))
        
        # Position nodes in layers
        for level, nodes in sorted(levels.items()):
            n_nodes = len(nodes)
            for i, node in enumerate(nodes):
                x = (i - n_nodes / 2) * 2.0  # Horizontal spread
                y = -level * 2.0  # Vertical drop per level
                pos[node] = (x, y)
        
        return pos

    
    def plot(
        self,
        figsize=(16, 12),
        node_color_attr='label',
        node_size_attr='count',
        outpath=None,
        colormap='tab20',
    ):
        """Plot network graph.
        
        Parameters:
            figsize: Figure size (width, height)
            node_color_attr: Node attribute for coloring ('label' uses cluster labels)
            node_size_attr: Node attribute for sizing ('count' uses collapse_counts)
            outpath: Path to save figure
            colormap: Matplotlib colormap name
        """
        if self.graph is None or self.pos is None:
            raise ValueError("Build graph and compute layout first")
        
        console.print("Plotting network...")
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Prepare node colors
        if node_color_attr == 'label' and self.labels is not None:
            node_colors = [self.labels[i] for i in self.graph.nodes()]
        else:
            node_colors = 'lightblue'
        
        # Prepare node sizes
        if node_size_attr == 'count':
            node_sizes = [
                300 + 100 * self.graph.nodes[i].get('count', 1)
                for i in self.graph.nodes()
            ]
        else:
            node_sizes = 300
        
        # Draw edges - distinguish MST edges from others
        mst_edges = [(u, v) for u, v in self.graph.edges() if self.graph[u][v].get('is_mst', False)]
        non_mst_edges = [(u, v) for u, v in self.graph.edges() if not self.graph[u][v].get('is_mst', False)]
        
        # Draw non-MST edges (lighter)
        if non_mst_edges:
            nx.draw_networkx_edges(
                self.graph,
                self.pos,
                edgelist=non_mst_edges,
                ax=ax,
                alpha=0.15,
                width=0.3,
                edge_color='gray',
                style='dotted'
            )
        
        # Draw MST edges (darker, thicker)
        if mst_edges:
            nx.draw_networkx_edges(
                self.graph,
                self.pos,
                edgelist=mst_edges,
                ax=ax,
                alpha=0.6,
                width=1.5,
                edge_color='darkblue'
            )
        
        # Draw nodes
        nx.draw_networkx_nodes(
            self.graph,
            self.pos,
            ax=ax,
            node_color=node_colors,
            node_size=node_sizes,
            cmap=colormap,
            alpha=0.8,
            edgecolors='black',
            linewidths=0.5,
        )
        
        # Draw labels (sparse to avoid clutter)
        labels_to_draw = {}
        n_nodes = self.graph.number_of_nodes()
        if n_nodes < 100:
            labels_to_draw = {i: self.graph.nodes[i].get('name', str(i)) for i in self.graph.nodes()}
        
        nx.draw_networkx_labels(
            self.graph,
            self.pos,
            labels_to_draw,
            ax=ax,
            font_size=6,
            font_color='black'
        )
        
        n_edges = self.graph.number_of_edges()
        ax.set_title(f"Network Graph ({n_nodes} nodes, {n_edges} edges)")
        ax.axis('off')
        plt.tight_layout()
        
        if outpath is not None:
            Path(outpath).parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(outpath, dpi=150, bbox_inches='tight')
            console.print(f"  [green]✓[/green] Plot saved to {outpath}")
        
        plt.close()


def run_network_visualization(
    distance_matrix,
    labels=None,
    sample_names=None,
    layout_method='spring',
    edge_threshold=None,
    edge_percentile=None,
    use_mst=True,
    mst_multiplier=1.5,
    outdir=None,
):
    """Complete network visualization workflow.
    
    Parameters:
        distance_matrix: Distance matrix
        labels: Cluster labels
        sample_names: Sample names
        collapse_counts: Sample count dict
        layout_method: 'spring', 'spectral', 'equidistant', 'circular'
        edge_threshold: Distance threshold for edges
        edge_percentile: Percentile-based edge threshold (0-100)
        use_mst: Use Minimum Spanning Tree for hierarchical structure (default: True)
        mst_multiplier: Multiplier for MST edges (1.0 = MST only, 1.5 = MST + additional context edges)
        outdir: Output directory
        
    Returns:
        dict with graph, positions, and output paths
    """
    net = NetworkGraph(
        distance_matrix,
        labels=labels,
        sample_names=sample_names,
    )
    
    net.build_graph(edge_threshold=edge_threshold, edge_percentile=edge_percentile)
    net.compute_layout(method=layout_method)
    
    outpath = None
    if outdir is not None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        outpath = str(outdir / f'network_{layout_method}.png')
    
    net.plot(outpath=outpath)
    
    return {
        'graph': net.graph,
        'positions': net.pos,
        'plot_path': outpath,
    }


__all__ = [
    'NetworkGraph',
    'run_network_visualization',
]

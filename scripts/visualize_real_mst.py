"""MST visualization for cgMLST data using spring or GrapeTree layout."""

import sys
import os
import tempfile
import math
from pathlib import Path
from typing import Dict, Iterable, Tuple, Set

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from rich.console import Console
import typer

# Keep your project import paths as before
sys.path.insert(0, str(Path(__file__).parent.parent))

from alleleatlas.core.input import detect_and_normalize_profile
from alleleatlas.core.distances import compute_distance_matrices

app = typer.Typer(help="MST visualization for cgMLST data")

console = Console()

# -----------------------
# Configuration / knobs
# -----------------------
DEFAULT_TEMP_NPROC = 4
DISTANCE_EDGE_CUTOFF = 20  # used for "short" edges when finding local components
DEFAULT_WINSORIZE_PCT = 99  # default winsorization percentile for grapetree layout
PLOT_OUTDIR = "out_real_mst"
FIGSIZE = (20, 16)
SEED = 42


# -----------------------
# Utilities
# -----------------------
def safe_remove(path: str):
    try:
        os.unlink(path)
    except Exception:
        pass


# -----------------------
# Data loading / distances
# -----------------------
def compute_distance_matrix(cgmlst_path: str, nproc: int = DEFAULT_TEMP_NPROC) -> np.ndarray:
    """
    Read cgMLST (PathogenWatch) profile, normalize to TSV via alleleatlas helper,
    call alleleatlas.compute_distance_matrices and return the distance matrix.
    """
    console.print("[bold]Loading and normalizing cgMLST data...[/bold]")
    df, _ = detect_and_normalize_profile(cgmlst_path, remap_alleles=True)
    console.print(f"  Profiles: {len(df)}, loci: {max(0, df.shape[1]-1)}")

    # write normalized TSV to temp file (alleleatlas expects file input)
    tf = tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False)
    normalized_path = tf.name
    try:
        df.to_csv(tf, sep="\t", index=False, header=True)
        tf.close()

        console.print("  Computing distance matrix (AlleleAtlas)...")
        distance_matrix, *_ = compute_distance_matrices(normalized_path, nproc=nproc)

    finally:
        safe_remove(normalized_path)

    # normalize possible 3D output => 2D symmetric matrix
    if distance_matrix.ndim == 3:
        # alleleatlas may return a lower-triangular block; attempt to convert to full symmetric
        arr = distance_matrix[:, :, 0].astype(float)
        n_samples = arr.shape[1]
        start = n_samples - arr.shape[0]
        sym = np.zeros((n_samples, n_samples), dtype=float)
        sym[start:, :] = arr
        for i in range(n_samples):
            for j in range(i):
                sym[j, i] = sym[i, j]
        distance_matrix = sym

    return distance_matrix


# -----------------------
# MST construction
# -----------------------
def build_mst_from_distances(distance_matrix: np.ndarray) -> nx.Graph:
    """Create a complete graph from distance matrix, then compute MST."""
    n = distance_matrix.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(n))

    # add undirected edges with 'weight' and 'distance' attrs
    for i in range(n):
        for j in range(i + 1, n):
            d = float(distance_matrix[i, j])
            if d >= 0 and not np.isnan(d):
                G.add_edge(i, j, weight=d, distance=d)

    mst = nx.minimum_spanning_tree(G, weight="weight")

    # basic node attributes
    for node in mst.nodes():
        mst.nodes[node]["name"] = f"S{node}"
        mst.nodes[node]["count"] = 1

    # mark MST edges
    for u, v in mst.edges():
        mst[u][v]["is_mst"] = True

    console.print(f"[bold]Built MST[/bold]: {mst.number_of_nodes()} nodes, {mst.number_of_edges()} edges")
    return mst


# -----------------------
# Component analysis
# -----------------------
def analyze_subgraph_metrics(G: nx.Graph, nodes: Iterable[int]) -> Dict:
    """Return simple metrics for a node set within G (density, internal edges, external connections)."""
    sub = G.subgraph(nodes)
    n = sub.number_of_nodes()
    if n <= 1:
        return {"size": n, "density": 1.0, "internal_edges": 0, "external_connections": 0}

    internal_edges = sub.number_of_edges()
    max_edges = n * (n - 1) / 2
    density = internal_edges / max_edges if max_edges > 0 else 0.0

    external = 0
    for node in sub.nodes():
        for nbr in G.neighbors(node):
            if nbr not in sub:
                external += 1

    return {"size": n, "density": density, "internal_edges": internal_edges, "external_connections": external}


# -----------------------
# GrapeTree layout (from te.py)
# -----------------------
def safe_asin(x: float) -> float:
    """Clamp to avoid numerical issues in asin."""
    return math.asin(max(-1.0, min(1.0, x)))


def polar_to_cartesian(r: float, theta: float) -> Tuple[float, float]:
    """Convert polar coords (r, theta in radians) to Cartesian (x, y)."""
    return (r * math.cos(theta), r * math.sin(theta))


class StaticGrapeTreeLayout:
    """
    Implements a static GrapeTree layout for a rooted tree (NetworkX MST).
    Each node is drawn as a circle with radius proportional to its weight.
    Edges become radial branches in a tree-like layout.
    """
    def __init__(self, G: nx.Graph, root: int = None, k: float = 1.0, default_branch_len: float = 1.0):
        self.G = G
        self.k = float(k)
        self.default_branch_len = default_branch_len
        self.root = root if root is not None else self._find_root()
        if self.root is None:
            self.root = list(G.nodes())[0] if G.number_of_nodes() > 0 else None
        
        # Build children dict using BFS from root
        self.children = {n: [] for n in G.nodes()}
        if self.root is not None:
            parent = {self.root: None}
            queue = [self.root]
            while queue:
                u = queue.pop(0)
                for v in G.neighbors(u):
                    if v not in parent:
                        parent[v] = u
                        self.children[u].append(v)
                        queue.append(v)
        
        self.r_node = {}
        self.pos = {}

    def _find_root(self) -> int:
        """Heuristic to find root node."""
        if self.G.number_of_nodes() == 0:
            return None
        # For undirected MST, just pick a node with highest degree or return first
        return max(self.G.nodes(), key=lambda n: self.G.degree(n))

    def node_radius(self, n: int) -> float:
        """Compute radius r based on node weight: r = weight**(k/2)"""
        if n in self.r_node:
            return self.r_node[n]
        m = self.G.nodes[n].get('weight', 1)
        r = float(m) ** (self.k / 2.0)
        self.r_node[n] = r
        return r

    def edge_length(self, parent: int, child: int) -> float:
        """Return branch length for edge parent->child."""
        if self.G.has_edge(parent, child):
            return self.G.edges[parent, child].get('length', self.default_branch_len)
        return self.default_branch_len

    def compute_dcs_recursive(self, v: int, s: float) -> Dict:
        """Compute DCS (Dynamic Central Span) recursively."""
        if len(self.children[v]) == 0:
            r = self.node_radius(v)
            return {'r': r, 'a': 0.0}

        r_ccs_list = []
        a_ccs_list = []
        for c in self.children[v]:
            dcs_c = self.compute_dcs_recursive(c, s)
            r_dcs_c, a_dcs_c = dcs_c['r'], dcs_c['a']
            l = self.edge_length(v, c)
            rc = self.node_radius(c)
            R_for_cs1 = max(1e-9, l + r_dcs_c)
            chord = max(1e-9, rc) * 2.0
            a1 = 2.0 * safe_asin(min(1.0, chord / (2.0 * R_for_cs1)))
            a2 = a_dcs_c
            a_ccs = max(a1, a2)
            r_ccs = l + r_dcs_c
            r_ccs_list.append(r_ccs)
            a_ccs_list.append(a_ccs)

        r_dcs = max(r_ccs_list) if r_ccs_list else self.node_radius(v)
        a_sum = sum(a_ccs_list)
        a_dcs = a_sum + s * len(a_ccs_list)
        return {'r': r_dcs, 'a': a_dcs}

    def find_s(self, max_trials: int = 30) -> float:
        """Iteratively search for optimal separating arc s."""
        V = list(self.G.nodes())
        s_t = 2.0 * math.pi / max(1, len(V))
        S_candidates = []
        for t in range(max_trials):
            max_a = 0.0
            for n in V:
                dcs = self.compute_dcs_recursive(n, s_t)
                max_a = max(max_a, dcs['a'])
            A_t = max_a
            if A_t <= 2.0 * math.pi + 1e-12:
                S_candidates.append((s_t, A_t))
            if A_t == 0:
                break
            s_t = s_t * (2.0 * math.pi / A_t)
        
        if not S_candidates:
            return 2.0 * math.pi / max(1, len(V)) * 0.1
        best = max(S_candidates, key=lambda x: x[1])
        return best[0]

    def place_nodes(self, s: float):
        """Compute final positions recursively."""
        def compute_all_dcs(node):
            dcs = {}
            def _rec(n):
                if len(self.children[n]) == 0:
                    d = {'r': self.node_radius(n), 'a': 0.0}
                    dcs[n] = d
                    return d
                r_ccs_list = []
                a_ccs_list = []
                ccs_items = []
                for c in self.children[n]:
                    d_child = _rec(c)
                    l = self.edge_length(n, c)
                    rc = self.node_radius(c)
                    R_for_cs1 = max(1e-9, l + d_child['r'])
                    chord = max(1e-9, rc) * 2.0
                    a1 = 2.0 * safe_asin(min(1.0, chord / (2.0 * R_for_cs1)))
                    a2 = d_child['a']
                    a_ccs = max(a1, a2)
                    r_ccs = l + d_child['r']
                    r_ccs_list.append(r_ccs)
                    a_ccs_list.append(a_ccs)
                    ccs_items.append({'child': c, 'a_ccs': a_ccs, 'r_ccs': r_ccs, 'l': l})
                r_dcs = max(r_ccs_list) if r_ccs_list else self.node_radius(n)
                a_dcs = sum(a_ccs_list) + s * len(a_ccs_list)
                d = {'r': r_dcs, 'a': a_dcs, 'ccs': ccs_items}
                dcs[n] = d
                return d
            _rec(node)
            return dcs

        if self.root is None:
            return
        
        dcs_all = compute_all_dcs(self.root)

        def _place(n, center_x, center_y, base_angle, angle_span):
            self.pos[n] = (center_x, center_y)
            info = dcs_all[n]
            if len(self.children[n]) == 0:
                return
            ccs = info.get('ccs', [])
            total = sum(item['a_ccs'] for item in ccs) + s * len(ccs)
            if total <= 0:
                total = max(1.0, float(len(ccs)))
                for i, item in enumerate(ccs):
                    child = item['child']
                    sub_angle = angle_span / max(1, len(ccs))
                    theta = base_angle + (i + 0.5) * sub_angle
                    l = item['l']
                    dx, dy = polar_to_cartesian(l, theta)
                    _place(child, center_x + dx, center_y + dy, theta - sub_angle / 2.0, sub_angle)
                return

            angle_cursor = base_angle
            for item in ccs:
                child = item['child']
                alloc = (item['a_ccs'] + s) / total * angle_span
                theta = angle_cursor + alloc / 2.0
                l = item['l']
                dx, dy = polar_to_cartesian(l, theta)
                _place(child, center_x + dx, center_y + dy, angle_cursor, alloc)
                angle_cursor += alloc

        root_dcs = dcs_all[self.root]
        full_span = root_dcs['a'] if root_dcs['a'] > 1e-12 else 2.0 * math.pi
        _place(self.root, 0.0, 0.0, -full_span / 2.0, full_span)

    def run(self, max_trials: int = 30) -> Dict[int, Tuple[float, float]]:
        """Compute layout and return pos dict node -> (x, y)"""
        s = self.find_s(max_trials=max_trials)
        self.place_nodes(s)
        # Apply gentle collision detection to reduce overlaps
        console.print("  Applying gentle collision detection...")
        return self.pos



# -----------------------
# Plotting
# -----------------------
def plot_graph(
    G: nx.Graph,
    positions: Dict[int, Tuple[float, float]],
    node_to_component: Dict[int, int],
    component_sizes: Dict[int, int],
    layout_name: str,
    distance_cutoff: float = DISTANCE_EDGE_CUTOFF,
    outdir: str = PLOT_OUTDIR,
    dataset_name: str = "dataset"
):
    """Plot graph using provided positions and save PNG to disk."""
    console.print(f"[bold]Plotting:[/bold] layout={layout_name}, nodes={G.number_of_nodes()}")

    fig, ax = plt.subplots(figsize=FIGSIZE, facecolor="white")
    ax.set_facecolor("white")

    # filter edges: hide very long edges or zero-length edges
    edges_to_draw = []
    for u, v in G.edges():
        d = G[u][v].get("weight", 0)
        if distance_cutoff > 0 and d > distance_cutoff:
            continue
        if d == 0:
            continue
        edges_to_draw.append((u, v))

    visible = G.edge_subgraph(edges_to_draw).copy()

    # degree for coloring
    degrees = dict(visible.degree()) if visible.number_of_nodes() > 0 else {}
    max_deg = max(degrees.values()) if degrees else 1

    # edge drawing with layout-specific styling
    if layout_name == "grapetree":
        # For GrapeTree, draw edges with varying alpha based on distance
        for u, v in edges_to_draw:
            d = G[u][v].get("weight", 0)
            # Closer edges (shorter distances) are more opaque
            alpha = 0.7 - (d / max(distance_cutoff, 1)) * 0.5 if distance_cutoff > 0 else 0.5
            alpha = max(0.1, min(0.7, alpha))
            nx.draw_networkx_edges(visible, positions, edgelist=[(u, v)], ax=ax, 
                                   alpha=alpha, width=0.8, edge_color="#555555")
    else:
        # Spring layout: draw all edges uniformly
        nx.draw_networkx_edges(visible, positions, ax=ax, alpha=0.3, width=0.6, edge_color="#222222")

    # node appearance
    node_colors = []
    node_sizes_px = []
    for n in G.nodes():
        deg = degrees.get(n, 0)
        if deg >= max_deg * 0.7:
            node_colors.append("#1E88E5")
        elif deg >= max_deg * 0.4:
            node_colors.append("#FFB300")
        else:
            node_colors.append("#424242")

        comp_id = node_to_component.get(n, None)
        comp_size = component_sizes.get(comp_id, 1) if comp_id is not None else 1
        # log scale for visual stability
        node_sizes_px.append(30 + np.log1p(comp_size) * 40)

    # glow layer
    nx.draw_networkx_nodes(G, positions, nodelist=list(G.nodes()), node_size=[s * 1.25 for s in node_sizes_px],
                           node_color=node_colors, alpha=0.18, linewidths=0)

    # main nodes
    nx.draw_networkx_nodes(G, positions, nodelist=list(G.nodes()), node_size=node_sizes_px,
                           node_color=node_colors, alpha=0.95, edgecolors="white", linewidths=0.9)

    # title
    n_original = G.number_of_nodes()
    n_components = len(set(node_to_component.values()))
    n_edges_total = G.number_of_edges()
    if distance_cutoff > 0:
        title = f"{dataset_name} — {layout_name} | {n_original} profiles ({n_components} islands) | {len(edges_to_draw)}/{n_edges_total} edges shown"
    else:
        title = f"{dataset_name} — {layout_name} | {n_original} profiles ({n_components} islands) | {n_edges_total} edges"
    ax.set_title(title, fontsize=18, fontweight="bold", pad=18)

    ax.axis("off")
    plt.tight_layout()

    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)
    outpath = outdir_p / f"{dataset_name}_{layout_name}.png"
    plt.savefig(str(outpath), dpi=150, bbox_inches="tight", facecolor="white")
    console.print(f"[green]✓[/green] Saved: {outpath}")
    plt.close(fig)


# -----------------------
# Orchestration: main run()
# -----------------------
def run_visualization(cgmlst_filepath: str, outdir: str = PLOT_OUTDIR, layout: str = "spring", winsorize_pct: float = None):
    """Internal function to run the visualization pipeline.
    
    Args:
        cgmlst_filepath: Path to cgMLST data
        outdir: Output directory for PNG
        layout: "spring" or "grapetree"
        winsorize_pct: For grapetree, cap edge lengths at this percentile (e.g., 95 caps at 95th percentile).
                       None = no winsorization
    """
    console.print(f"\n[bold]Real Data MST Visualization ({layout} layout)[/bold]\n")

    dataset_name = Path(cgmlst_filepath).stem
    dataset_name = f"mst_{dataset_name}"

    # 1) distances
    dist_matrix = compute_distance_matrix(cgmlst_filepath)

    # 2) MST
    G = build_mst_from_distances(dist_matrix)

    # 3) detect short-edge components (local islands)
    short_cutoff = DISTANCE_EDGE_CUTOFF
    short_edges = [(u, v) for u, v in G.edges() if float(G[u][v].get("weight", 0)) < short_cutoff]
    short_sub = G.edge_subgraph(short_edges)
    components_list = list(nx.connected_components(short_sub))

    # map nodes -> component id; nodes not in any short component become their own singletons
    node_to_component: Dict[int, int] = {}
    components: Dict[int, Set[int]] = {}
    for cid, comp_nodes in enumerate(components_list):
        components[cid] = set(comp_nodes)
        for n in comp_nodes:
            node_to_component[n] = cid

    # remaining nodes as singleton components
    for n in G.nodes():
        if n not in node_to_component:
            cid = len(components)
            components[cid] = {n}
            node_to_component[n] = cid

    # analyze & display simple stats
    component_sizes = {cid: len(nodes) for cid, nodes in components.items()}
    console.print(f"  Components detected: {len(components)} (short-edge cutoff={short_cutoff})")
    # print a few components summary
    top = sorted(component_sizes.items(), key=lambda x: x[1], reverse=True)[:6]
    for cid, size in top:
        m = analyze_subgraph_metrics(G, components[cid])
        console.print(f"    comp {cid}: size={size}, density={m['density']:.2f}, ext_conn={m['external_connections']}")

    # 4) compute layout
    if layout == "grapetree":
        console.print("[bold]Computing GrapeTree layout...[/bold]")
        # Set node weights for GrapeTree (component size)
        for n in G.nodes():
            comp_id = node_to_component.get(n, None)
            G.nodes[n]['weight'] = component_sizes.get(comp_id, 1)
        
        # Set edge lengths based on MST weights with optional winsorization
        edge_lengths = []
        for u, v in G.edges():
            edge_lengths.append(G[u][v].get('weight', 1.0))
        
        if winsorize_pct is not None and 0 < winsorize_pct < 100:
            # Cap edge lengths at the specified percentile
            max_len = np.percentile(edge_lengths, winsorize_pct)
            console.print(f"  Winsorizing edge lengths at {winsorize_pct}th percentile (cap={max_len:.2f})")
            for u, v in G.edges():
                orig_len = G[u][v].get('weight', 1.0)
                capped_len = min(orig_len, max_len)
                G[u][v]['length'] = capped_len
        else:
            for u, v in G.edges():
                G[u][v]['length'] = G[u][v].get('weight', 1.0)
        
        layouter = StaticGrapeTreeLayout(G, k=1.0, default_branch_len=1.0)
        pos = layouter.run(max_trials=30)
    else:  # spring (default)
        console.print("[bold]Computing spring layout...[/bold]")
        pos = nx.spring_layout(G, seed=SEED, weight="weight", k=4.5, iterations=4000)
    
    plot_graph(G, pos, node_to_component, component_sizes, layout_name=layout,
               distance_cutoff=DISTANCE_EDGE_CUTOFF, outdir=outdir, dataset_name=dataset_name)

    console.print("\n[bold green]✓ Complete![/bold green]")


@app.command()
def visualize(
    cgmlst_path: str = typer.Argument(
        ...,
        help="Path to cgMLST profile file (CSV, TSV, or gzipped)"
    ),
    output_dir: str = typer.Option(
        PLOT_OUTDIR,
        "--output-dir", "-o",
        help="Directory to save output PNG files"
    ),
    edge_cutoff: int = typer.Option(
        DISTANCE_EDGE_CUTOFF,
        "--edge-cutoff", "-e",
        help="Distance cutoff for detecting local components (short edges)"
    ),
    layout: str = typer.Option(
        "grapetree",
        "--layout", "-l",
        help="Layout algorithm: 'spring' (force-directed) or 'grapetree' (tree-like circular)"
    ),
    winsorize: float = typer.Option(
        DEFAULT_WINSORIZE_PCT,
        "--winsorize", "-w",
        help="For grapetree layout: cap edge lengths at this percentile (0-100). E.g., 95 caps at 95th percentile. Default is 99"
    ),
    seed: int = typer.Option(
        SEED,
        "--seed", "-s",
        help="Random seed for layout reproducibility"
    ),
) -> None:
    """Visualize cgMLST data as a minimum spanning tree.
    
    Generates spring layout (Fruchterman-Reingold force-directed) or GrapeTree visualization.
    """
    # Override config with CLI args
    global DISTANCE_EDGE_CUTOFF, SEED, PLOT_OUTDIR
    DISTANCE_EDGE_CUTOFF = edge_cutoff
    SEED = seed
    PLOT_OUTDIR = output_dir
    
    if layout not in ("spring", "grapetree"):
        console.print(f"[bold red]✗ Error:[/bold red] Invalid layout '{layout}'. Choose 'spring' or 'grapetree'")
        raise typer.Exit(code=1)
    
    try:
        run_visualization(cgmlst_path, outdir=output_dir, layout=layout, winsorize_pct=winsorize)
    except FileNotFoundError:
        console.print(f"[bold red]✗ Error:[/bold red] File not found: {cgmlst_path}")
        raise typer.Exit(code=1)
    except Exception as e:
        console.print(f"[bold red]✗ Error:[/bold red] {str(e)}")
        raise typer.Exit(code=1)


# -----------------------
# CLI entrypoint
# -----------------------
if __name__ == "__main__":
    app()

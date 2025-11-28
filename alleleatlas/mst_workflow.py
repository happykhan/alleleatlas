"""Simple MST workflow that mirrors the legacy getDistances/MSTrees behaviour.

Features:
- load and clean profile files (simple subset of legacy parser)
- compute distance matrices (symmetric, asymmetric, blockwise) using same formulas
- compute MST using networkx (symmetric -> undirected MST, asymmetric -> arborescence fallback)
- compute distance distribution, cluster by threshold and collapse nodes into supernodes
- plot MST with spring layout and optional convex hulls for clusters
"""
from __future__ import annotations
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib import colors as mcolors
from matplotlib import colors as mcolors
import os
import re
import numba as nb

try:
    from scipy.spatial import ConvexHull
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False


MISSING_TOKENS = {'0', 'N', '-'}


def read_profile(path: str):
    """Parse a profile file using the same simple logic as `MSTrees.backend`.

    Returns (names, profiles) where `profiles` is an ndarray of strings
    (alleles) and `names` is a list of sample names / STs.
    """
    names, profiles = [], []
    try:
        if path[-3:].lower().endswith('.gz'):
            fin = open(path, 'rt').readlines() if os.path.isfile(path) else path.split('\n')
        else:
            fin = open(path).readlines() if os.path.isfile(path) else path.split('\n')
    except Exception:
        fin = path.split('\n')

    allele_cols = None
    fmt = 'profile'
    for line_id, line in enumerate(fin):
        if line.startswith('#'):
            if not line.startswith('##'):
                header = line.strip().split('\t')
                allele_cols = [id for id, col in enumerate(header) if id > 0 and not col.startswith('#') and col.lower() not in {'st_id', 'st'}]
            continue
        if line.startswith('>'):
            fmt = 'fasta'
        else:
            fmt = 'profile'
            if allele_cols is None:
                header = line.strip().split('\t')
                allele_cols = [id for id, col in enumerate(header) if id > 0 and not col.startswith('#') and col.lower() not in {'st_id', 'st'}]
                line_id += 1
        break

    if fmt == 'fasta':
        for line in fin[line_id:]:
            if line.startswith('>'):
                names.append(line[1:].strip().split()[0])
                profiles.append([])
            else:
                profiles[-1].extend(line.strip().split())
        for id, p in enumerate(profiles):
            profiles[id] = list(''.join(p))
    else:
        for line in fin[line_id:]:
            part = line.strip().split('\t')
            if not part or not part[0]:
                continue
            names.append(part[0])
            if allele_cols is not None:
                profiles.append([part[i] if i < len(part) else '0' for i in allele_cols])
            else:
                profiles.append(part[1:])

    profiles = np.char.upper(np.array(profiles, dtype=str))
    names = [re.sub(r'[\(\)\ \,\"\';]', '_', n) for n in names]
    return np.array(names), profiles


def nonredundant(names, profiles):
    """Ported from old_code/MSTrees.nonredundant

    Expects names: np.array of strings, profiles: 2D array-like of tokens
    (strings or ints). Returns (names, profiles, embeded) where profiles
    are integer-encoded with missing -> 0 and duplicates collapsed into
    embeded groups.
    """
    # Ensure numpy arrays of strings for encoding logic similar to legacy
    names = np.array(names)
    # If profiles are numeric already, convert to strings for unique mapping
    prof_str = np.array(profiles, dtype=str)
    encoded_profile = np.array([np.unique(p, return_inverse=True)[1] + 1 for p in prof_str.T]).T
    # legacy treats '0'|'N'|'-' as missing -> set to 0
    mask_missing = (prof_str == '0') | (prof_str == 'N') | (prof_str == '-')
    encoded_profile[mask_missing] = 0

    # TODO: params['handle_missing'] isn't available here; we won't perform complete_delete
    names = names[np.lexsort(encoded_profile.T)]
    profiles_enc = encoded_profile[np.lexsort(encoded_profile.T)]
    presence = (np.sum(profiles_enc > 0, 1) > 0)
    names = names[presence]
    profiles_enc = profiles_enc[presence]

    uniqueness = np.concatenate([[1], np.sum(np.diff(profiles_enc, axis=0) != 0, 1) > 0])

    embeded = {names[0]: []}
    embeded_group = embeded[names[0]]
    for n, u in zip(names, uniqueness):
        if u == 0:
            embeded_group.append(n)
        else:
            embeded[n] = [n]
            embeded_group = embeded[n]
    names = names[uniqueness > 0]
    profiles_enc = profiles_enc[uniqueness > 0]
    return names, profiles_enc.astype(np.int32), embeded


def distance_symmetric(profiles: np.ndarray, handle_missing='pair_delete') -> np.ndarray:
    n_profiles, n_loci = profiles.shape
    if handle_missing in ('as_allele',):
        presences = np.ones_like(profiles, dtype=bool)
    elif handle_missing in ('pair_delete', 'absolute_distance'):
        presences = profiles > 0
    else:
        presences = np.repeat(np.sum(profiles > 0, 0) >= profiles.shape[0], profiles.shape[0]).reshape([profiles.shape[1], profiles.shape[0]]).T

    dist = np.zeros((n_profiles, n_profiles), dtype=np.float32)

    if handle_missing in ('pair_delete',):
        for i in range(n_profiles):
            profile = profiles[i]
            presence = presences[i]
            if i > 0:
                comparable = presences[:i] & presence
                al = np.sum(comparable, axis=1)
                diffs = np.sum((profiles[:i] != profile) & comparable, axis=1)
                adj = (diffs + 0.01) * float(presence.size) / (al + 0.01)
                dist[:i, i] = adj
                dist[i, :i] = adj
    else:
        for i in range(n_profiles):
            profile = profiles[i]
            presence = presences[i]
            if i > 0:
                comparable = presences[:i] & presence
                diffs = np.sum((profiles[:i] != profile) & comparable, axis=1)
                dist[:i, i] = diffs
                dist[i, :i] = diffs
    return dist


def distance_asymmetric(profiles: np.ndarray, handle_missing='pair_delete') -> np.ndarray:
    # legacy dual_dist expects an int matrix with 0 meaning missing
    # Two accepted input shapes:
    # - integer-encoded matrix where column 0 is ST/dummy and allele data starts at col 1
    #   In that case we call dual_dist on mat[:,1:] to match legacy mat[:,1:] semantics
    # - string token matrix (profiles) where we must map tokens to per-column integers
    n_profiles, n_loci = profiles.shape

    # Case A: already integer-encoded matrix
    if np.issubdtype(profiles.dtype, np.integer):
        mat = profiles
        # if there is a leading dummy/ST column emulate legacy by slicing it off
        if mat.shape[1] >= 2:
            allele_mat = mat[:, 1:]
        else:
            allele_mat = mat.copy()
        try:
            d3 = dual_dist(allele_mat, 0, n_profiles, 0.05)
            dist = np.zeros((n_profiles, n_profiles), dtype=np.float32)
            for i in range(n_profiles):
                for j in range(i):
                    dist[i, j] = d3[i, j, 0]
                    dist[j, i] = d3[i, j, 1]
            return dist
        except Exception:
            # fall through to python fallback on the integer matrix
            presences = mat > 0
            dist = np.zeros((n_profiles, n_profiles), dtype=np.float32)
            if handle_missing not in ('absolute_distance',):
                for i in range(n_profiles):
                    profile = mat[i]
                    presence = presences[i]
                    diffs = np.sum(((mat != profile) & presence), axis=1) * float(presence.size) / (np.sum(presence) + 1e-12)
                    dist[:, i] = diffs
            else:
                for i in range(n_profiles):
                    profile = mat[i]
                    presence = presences[i]
                    diffs = np.sum((mat != profile) & presence, axis=1)
                    dist[:, i] = diffs
            return dist

    # Case B: string token matrix -> build integer matrix per-column (0 means missing)
    mat = np.zeros((n_profiles, n_loci), dtype=np.int32)
    for c in range(n_loci):
        col = profiles[:, c].astype(str)
        unique, inv = np.unique(col, return_inverse=True)
        # treat '0','N','-' as missing (0)
        map_vals = np.arange(1, unique.size + 1, dtype=np.int32)
        for i_u, u in enumerate(unique):
            if u in ('0', 'N', '-'):
                map_vals[i_u] = 0
        mat[:, c] = map_vals[inv]

    try:
        # call ported dual_dist on the allele-only integer matrix
        d3 = dual_dist(mat, 0, n_profiles, 0.05)
        dist = np.zeros((n_profiles, n_profiles), dtype=np.float32)
        for i in range(n_profiles):
            for j in range(i):
                dist[i, j] = d3[i, j, 0]
                dist[j, i] = d3[i, j, 1]
        # Legacy `methods.distance` divides by number of loci for non-absolute distances
        if handle_missing != 'absolute_distance':
            n_loci = mat.shape[1]
            return dist.astype(np.float32) / float(n_loci)
        return dist
    except Exception:
        # fallback to simple python implementation
        presences = mat > 0
        dist = np.zeros((n_profiles, n_profiles), dtype=np.float32)
        if handle_missing not in ('absolute_distance',):
            for i in range(n_profiles):
                profile = mat[i]
                presence = presences[i]
                diffs = np.sum(((mat != profile) & presence), axis=1) * float(presence.size) / (np.sum(presence) + 1e-12)
                dist[:, i] = diffs
        else:
            for i in range(n_profiles):
                profile = mat[i]
                presence = presences[i]
                diffs = np.sum((mat != profile) & presence, axis=1)
                dist[:, i] = diffs
        return dist


def distance_blockwise(profiles: np.ndarray, handle_missing=0.01) -> np.ndarray:
    if isinstance(handle_missing, str):
        handle_missing = float(handle_missing)
    presences = profiles > 0
    n_profiles = profiles.shape[0]
    distances = np.zeros((n_profiles, n_profiles), dtype=np.float32)
    for i in range(n_profiles):
        profile = profiles[i]
        diffs = np.hstack([np.zeros([n_profiles, 1], dtype=int), profiles - profile, np.zeros([n_profiles, 1], dtype=int)])
        d1 = np.sum((diffs[:, 1:] != diffs[:, :-1]) & (diffs[:, 1:] != 0), 1)
        d2 = np.sum(diffs != 0, 1) - d1
        distances[:, i] = (d1 + d2 * handle_missing)
    return distances


@nb.jit(nopython=True)
def dual_dist(mat, s, e, allowed_missing=0.05):
    dist = np.zeros((e - s, mat.shape[0], 2), dtype=np.int32)
    n_loci = mat.shape[1]
    for i in range(s, e):
        ql = 0
        for kk in range(n_loci):
            if mat[i, kk] > 0:
                ql += 1
        for j in range(i):
            rl = 0.0
            ad = 1e-4
            al = 1e-4
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
    dist = np.zeros((e-s, mat.shape[0], 2), dtype=np.int32 )
    n_loci = mat.shape[1]
    for i in range(s, e) :
        for j in range(i) :
            ad, al = 0., 0.
            for k in range(n_loci) :
                if mat[j, k] > 0 :
                    if mat[i, k] > 0 :
                        al += 1
                        if mat[i, k] != mat[j, k] :
                            ad += 1
            dist[i-s, j, 0] = int( -np.log(1.-(ad+0.5)/(al+1.0)) * n_loci * 100. + 0.5)
    return dist


def heuristic_weights(dist: np.ndarray, names_embedded: list, method: str = 'harmonic') -> np.ndarray:
    # names_embedded is used to compute n_str (size of groups)
    n_str = np.array([len(g) for g in names_embedded])
    if method == 'harmonic':
        weights = dist.shape[0] / np.sum(1.0 / (dist + 0.1), axis=1)
        cw = np.vstack([-np.array(n_str), weights])
        weights[np.lexsort(cw)] = np.arange(dist.shape[0], dtype=float) / dist.shape[0]
        return weights
    else:  # eBurst
        weights = np.apply_along_axis(lambda r: np.bincount(r.astype(int), minlength=int(np.max(dist).astype(int) + 2)), 1, np.hstack([dist.astype(int), np.array([[int(np.max(dist)) + 1]] * dist.shape[1])]))
        weights = weights.T
        weights[0] += n_str
        dist_order = np.concatenate([[0], np.arange(weights.shape[1] - 1, 0, -1)])
        orders = np.lexsort(-weights.T[dist_order])
        w = np.zeros(dist.shape[0])
        w[orders] = (np.arange(orders.size)) / float(orders.size)
        return w


def build_mst(dist: np.ndarray, weights: np.ndarray, matrix_type='asymmetric') -> list:
    # Mirror methods._symmetric/_asymmetric in behaviour (simplified)
    if matrix_type == 'blockwise':
        # scaled variant in legacy code
        scaled = np.round(dist * 10000.0, 0) + weights.reshape([weights.size, -1])
        g = nx.Graph()
        n = scaled.shape[0]
        for i in range(n):
            for j in range(i):
                g.add_edge(i, j, weight=int(scaled[i, j]))
        ms = nx.minimum_spanning_tree(g)
        return [[u, v, int(d['weight']) / 10000.0] for u, v, d in ms.edges(data=True)]
    elif matrix_type == 'symmetric':
        d2 = np.round(dist, 0) + weights.reshape([weights.size, -1])
        np.fill_diagonal(d2, 0.0)
        d2[d2 > d2.T] = d2.T[d2 > d2.T]
        g = nx.Graph()
        n = d2.shape[0]
        for i in range(n):
            for j in range(i):
                g.add_edge(i, j, weight=int(d2[i, j]))
        ms = nx.minimum_spanning_tree(g)
        return [[u, v, int(d['weight'])] for u, v, d in ms.edges(data=True)]
    else:  # asymmetric
        d2 = np.round(dist, 0) + weights.reshape([weights.size, -1])
        np.fill_diagonal(d2, 0.0)
        g = nx.DiGraph()
        n = d2.shape[0]
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                g.add_edge(i, j, weight=int(d2[i, j]))
        try:
            ms = nx.minimum_spanning_arborescence(g)
            return [[u, v, int(d['weight'])] for u, v, d in ms.edges(data=True)]
        except Exception:
            # fallback: undirected MST on symmetric projection
            gu = nx.Graph()
            for i in range(n):
                for j in range(i):
                    w = min(int(d2[i, j]), int(d2[j, i]))
                    gu.add_edge(i, j, weight=w)
            ms = nx.minimum_spanning_tree(gu)
            return [[u, v, int(d['weight'])] for u, v, d in ms.edges(data=True)]


def collapse_components(edges: list, cutoff: float):
    """Collapse nodes connected by edges with weight < cutoff into super-nodes.

    Returns mapping node->component_id and a new list of edges between components with min weight.
    """
    g = nx.Graph()
    for u, v, w in edges:
        g.add_edge(u, v, weight=w)
    merge_g = nx.Graph(((u, v) for u, v, d in g.edges(data=True) if d['weight'] < cutoff))
    comps = list(nx.connected_components(merge_g))
    node2comp = {}
    for cid, comp in enumerate(comps):
        for n in comp:
            node2comp[n] = cid
    # nodes that were isolated
    all_nodes = set(g.nodes())
    for n in all_nodes:
        if n not in node2comp:
            node2comp[n] = len(comps)
            comps.append({n})

    # build inter-component edges with min weight
    comp_edges = {}
    for u, v, w in edges:
        cu, cv = node2comp[u], node2comp[v]
        if cu == cv:
            continue
        key = tuple(sorted((cu, cv)))
        comp_edges[key] = min(comp_edges.get(key, float('inf')), w)

    comp_edge_list = [[a, b, w] for (a, b), w in comp_edges.items()]
    return node2comp, comp_edge_list, comps


def plot_mst(edges: list, names: list, node2comp=None, comps=None, layout='spring', iterations=200, figsize=(10, 8), out_file='mst.png'):
    # build graph
    G = nx.Graph()
    for u, v, w in edges:
        G.add_edge(u, v, weight=w)
    names = list(names)
    if layout == 'spring':
        pos = nx.spring_layout(G, iterations=iterations, seed=42)
    elif layout == 'grapetree':
        # convert NetworkX edges to list of (u,v,weight)
        eds = []
        for u, v, data in G.edges(data=True):
            w = data.get('weight', 1.0)
            eds.append((u, v, w))
        pos = grapetree_layout(eds, root=0)
    else:
        pos = nx.spring_layout(G, iterations=iterations, seed=42)

    plt.figure(figsize=figsize)
    # node colors by component
    if node2comp is None:
        colors = [mcolors.to_hex(plt.get_cmap('tab20')(0)) for _ in G.nodes()]
    else:
        comp_ids = [node2comp.get(n, -1) for n in G.nodes()]
        unique = sorted(set(comp_ids))
        cmap = plt.get_cmap('tab20')
        color_map = {cid: mcolors.to_hex(cmap(i % 20)) for i, cid in enumerate(unique)}
        colors = [color_map[cid] for cid in comp_ids]

    weights = [d['weight'] for _, _, d in G.edges(data=True)]
    nx.draw_networkx_edges(G, pos, alpha=0.6)
    nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=80)
    nx.draw_networkx_labels(G, pos, {i: names[i] for i in G.nodes()}, font_size=6)

    # convex hulls for comps if available
    if _HAVE_SCIPY and comps is not None:
        for cid, comp in enumerate(comps):
            if len(comp) < 3:
                continue
            pts = np.array([pos[n] for n in comp])
            try:
                hull = ConvexHull(pts)
                poly = pts[hull.vertices]
                plt.fill(poly[:, 0], poly[:, 1], alpha=0.1)
            except Exception:
                pass

    plt.axis('off')
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close()
    return out_file


def grapetree_layout(edges: list, root=0, node_size=None, link_value_scale=1.0):
    """A simplified port of the GrapeTree greedy radial layout.

    Produces a positions dict {node: (x,y)} for the MST edges. This is
    a lightweight approximation of the JS greedy_layout sufficient for
    static plotting.
    """
    # build adjacency
    adj = {}
    for u, v, w in edges:
        adj.setdefault(u, []).append((v, w))
        adj.setdefault(v, []).append((u, w))

    # build rooted tree using BFS
    parent = {root: None}
    children = {root: []}
    link_value = {}
    from collections import deque
    q = deque([root])
    while q:
        n = q.popleft()
        for (nbr, w) in adj.get(n, []):
            if nbr in parent:
                continue
            parent[nbr] = n
            children.setdefault(n, []).append(nbr)
            children.setdefault(nbr, [])
            link_value[(n, nbr)] = w
            q.append(nbr)

    # compute subtree sizes
    size = {}
    def compute_size(n):
        s = 1
        for c in children.get(n, []):
            s += compute_size(c)
        size[n] = s
        return s
    compute_size(root)

    # assign angular spans recursively
    positions = {}
    def assign(n, angle, span, radius=0.0):
        # place node at polar coords (radius, angle)
        positions[n] = (radius * np.cos(angle), radius * np.sin(angle))
        # distribute child spans
        total = sum(size[c] for c in children.get(n, []))
        if total == 0:
            return
        start = angle - span/2
        for c in children.get(n, []):
            frac = size[c] / total
            child_span = span * frac
            child_angle = start + child_span/2
            # radial increment proportional to link weight if available
            w = link_value.get((n, c), 1.0)
            assign(c, child_angle, child_span, radius + np.log(1 + w) * 0.15)
            start += child_span

    assign(root, 0.0, np.pi*2, radius=0.0)
    return positions

    plt.figure(figsize=figsize)
    # node colors by component
    if node2comp is None:
        colors = 'C0'
    else:
        comp_ids = [node2comp.get(n, -1) for n in G.nodes()]
        unique = sorted(set(comp_ids))
        cmap = plt.get_cmap('tab20')
        color_map = {cid: cmap(i % 20) for i, cid in enumerate(unique)}
        colors = [color_map[cid] for cid in comp_ids]

    weights = [d['weight'] for _, _, d in G.edges(data=True)]
    nx.draw_networkx_edges(G, pos, alpha=0.6)
    nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=80)
    nx.draw_networkx_labels(G, pos, {i: names[i] for i in G.nodes()}, font_size=6)

    # convex hulls for comps if available
    if _HAVE_SCIPY and comps is not None:
        for cid, comp in enumerate(comps):
            if len(comp) < 3:
                continue
            pts = np.array([pos[n] for n in comp])
            try:
                hull = ConvexHull(pts)
                poly = pts[hull.vertices]
                plt.fill(poly[:, 0], poly[:, 1], alpha=0.1)
            except Exception:
                pass

    plt.axis('off')
    plt.tight_layout()
    plt.savefig(out_file, dpi=150)
    plt.close()
    return out_file


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('profile', help='profile file to run')
    parser.add_argument('--matrix', choices=['symmetric', 'asymmetric', 'blockwise'], default='asymmetric')
    parser.add_argument('--heuristic', choices=['harmonic', 'eBurst'], default='harmonic')
    parser.add_argument('--collapse', type=float, default=1.0, help='edge weight cutoff to collapse into supernodes (edges < cutoff merged)')
    parser.add_argument('--plot', default='mst.png')
    args = parser.parse_args()

    names, profiles = read_profile(args.profile)

    # Apply legacy-style nonredundant preprocessing (collapses identical rows)
    names_arr, profiles_enc, embeded = nonredundant(names, profiles)

    # For distance computations we use the encoded profiles
    profiles_to_use = profiles_enc

    if args.matrix == 'symmetric':
        dist = distance_symmetric(profiles_to_use)
    elif args.matrix == 'blockwise':
        dist = distance_blockwise(profiles_to_use)
    else:
        # legacy dual_dist expects a matrix with first column being an ST column
        # and the allele data starting at column 1; emulate that by prefixing
        # a dummy column so mat[:,1:] yields the allele matrix
        dummy = np.zeros((profiles_to_use.shape[0], 1), dtype=profiles_to_use.dtype)
        mat_for_dual = np.hstack([dummy, profiles_to_use])
        # call the ported dual_dist directly via distance_asymmetric wrapper
        dist = distance_asymmetric(mat_for_dual)

    # embedded groups: every name individual (no collapse by identical rows here)
    embeded = [[n] for n in names]
    weights = heuristic_weights(dist, embeded, method=args.heuristic)
    edges = build_mst(dist, weights, matrix_type=args.matrix)

    node2comp, comp_edges, comps = collapse_components(edges, args.collapse)
    out = plot_mst(edges, names, node2comp=node2comp, comps=comps, out_file=args.plot)
    print('Wrote', out)

"""
Enhanced UMAP/MDS visualization that preserves sample counts from collapsed profiles.

When working with collapsed profiles (where each cluster represents multiple original
samples), this module:
- Reads sample counts from collapse metadata
- Scales node sizes by sample count in visualizations
- Propagates counts to downstream MST/network visualizations
"""

from pathlib import Path
import json
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from alleleatlas.cluster.pHierCC import prepare_mat
from alleleatlas.cluster.getDistance import getDistance
from alleleatlas.select_hc_breakpoints import parse_eval_tsv, select_breakpoints
from alleleatlas.visualization.utils import num_from, norm_id


def _load_collapse_metadata(profile_path):
    """
    Try to load collapse metadata for a profile.
    
    Searches for metadata file in same directory with pattern:
    collapse_metadata_*_*.json
    
    Returns dict or None if not found.
    """
    profile_dir = Path(profile_path).parent
    
    # Try exact match in same directory
    for meta_file in profile_dir.glob('collapse_metadata_*.json'):
        try:
            with open(meta_file) as f:
                return json.load(f)
        except Exception:
            pass
    
    return None


def plot_emb(emb, out_path, title, names, clabels, top, colors_list, plotted_colors, hc_col, sample_counts=None):
    """
    Plot embedding (UMAP or MDS).
    
    If sample_counts provided, scale node sizes by count.
    """
    figsize = (10, 8) if sample_counts else (8, 6)
    plt.figure(figsize=figsize)
    
    # Scale node size by sample count if provided
    if sample_counts is not None:
        # Get count for each sample (assume names match cluster IDs or indices)
        sizes = []
        for i, name in enumerate(names):
            # Try to get count by cluster ID (for collapsed profiles)
            count = sample_counts.get(str(i), 1)
            if count is None:
                count = 1
            sizes.append(max(20, min(300, count * 0.5)))  # Scale between 20-300
        s = sizes
    else:
        s = 24
    
    plt.scatter(emb[:, 0], emb[:, 1], c=plotted_colors, s=s, alpha=0.7, edgecolors='k', linewidth=0.5)
    
    for i_lab, lab in enumerate(top):
        idxs = [i for i, v in enumerate(clabels) if v == lab]
        if not idxs:
            continue
        cx = emb[idxs, 0].mean()
        cy = emb[idxs, 1].mean()
        plt.scatter([cx], [cy], s=200, marker='X', color=colors_list[i_lab], edgecolors='k', linewidth=1.5, zorder=10)
    
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker='o', color='w', label=str(lab), markerfacecolor=colors_list[i], markersize=8, markeredgecolor='k') for i, lab in enumerate(top)]
    other_handle = Line2D([0], [0], marker='o', color='w', label='Other', markerfacecolor='#d0d0d0', markersize=8, markeredgecolor='k')
    handles.append(other_handle)
    
    if sample_counts is not None:
        # Add size legend
        size_legend_text = "Node size ∝ sample count"
        plt.text(0.02, 0.98, size_legend_text, transform=plt.gca().transAxes,
                fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    if handles:
        plt.legend(handles=handles, title=f'{hc_col} (top groups)', loc='upper right')
    
    plt.title(title, fontsize=12, fontweight='bold')
    plt.xlabel(title + '1')
    plt.ylabel(title + '2')
    plt.tight_layout()
    plt.savefig(str(out_path), dpi=150, bbox_inches='tight')
    plt.close()


def run(profile, outdir, nproc=2, metadata=None):
    """
    Run UMAP/MDS with optional collapse metadata for sample counts.
    
    Parameters:
        profile (str): Path to profile file
        outdir (str): Output directory
        nproc (int): Number of processes
        metadata (dict): Optional collapse metadata with sample_counts. If None, tries to auto-load.
    """
    PROFILE = Path(profile)
    OUT = Path(outdir)
    OUT.mkdir(parents=True, exist_ok=True)

    mat, names = prepare_mat(str(PROFILE))

    # Load collapse metadata if available
    if metadata is None:
        metadata = _load_collapse_metadata(str(PROFILE))
    
    sample_counts = None
    if metadata is not None and 'sample_counts' in metadata:
        sample_counts = metadata['sample_counts']
        print(f"Loaded collapse metadata: {len(sample_counts)} clusters with sample counts")

    # compute distances (for MDS)
    from multiprocessing import Pool
    pool = Pool(nproc)
    try:
        dist3 = getDistance(mat, 'dual_dist', pool, start=0, allowed_missing=0.03)
    finally:
        pool.close()
        pool.join()

    D = dist3[:, :, 0].astype(float)
    D = D + D.T
    np.fill_diagonal(D, 0.0)

    # UMAP from features
    embedding_umap = None
    try:
        import umap
        features = mat[:, 1:]
        reducer = umap.UMAP(metric='euclidean', n_components=2, n_neighbors=15, min_dist=0.1, random_state=42)
        embedding_umap = reducer.fit_transform(features)
    except Exception as e:
        print(f"UMAP failed: {e}")
        embedding_umap = None

    # classical MDS
    mds_emb = None
    try:
        from sklearn.manifold import MDS
        mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42, n_init=4, max_iter=500)
        mds_emb = mds.fit_transform(D)
    except Exception as e:
        print(f"MDS failed: {e}")
        mds_emb = None

    # choose HC column for coloring (try common eval file locations)
    sil = sim = labels = None
    for cand in [OUT / 'hierCC.eval', OUT / 'hierCC.eval.tsv', Path('output/hierCC.eval')]:
        if cand.exists():
            sil, sim, labels = parse_eval_tsv(str(cand))
            break

    chosen = []
    if sil is not None:
        _, _, _, chosen = select_breakpoints(sil, sim, top_n=10)

    best = None
    best_val = None
    for hc_level, idx, score, s, n in chosen:
        lbl = labels[idx] if (labels is not None and idx < len(labels)) else None
        num = num_from(lbl)
        val = num if num is not None else idx
        if best_val is None or val < best_val:
            best_val = val
            best = (hc_level, idx, lbl)

    # read hierCC table
    hc_col = None
    hdf = None
    hpath = None
    for cand in [OUT / 'hierCC.HierCC.gz', OUT / 'hierCC.HierCC', Path('output/hierCC.HierCC.gz')]:
        if cand.exists():
            hpath = cand
            break

    if hpath is not None:
        compression = 'gzip' if hpath.suffix == '.gz' else None
        hdf = pd.read_csv(hpath, sep='\t', compression=compression, dtype=str)
        sample_col = hdf.columns[0]
        hc_cols = [c for c in hdf.columns if str(c).upper().startswith('HC')]
        if best is not None:
            _, idx, lbl = best
            if lbl and lbl in hdf.columns:
                hc_col = lbl
            else:
                tgt = num_from(lbl) if lbl else None
                if tgt is not None:
                    for c in hc_cols:
                        if num_from(c) == tgt:
                            hc_col = c
                            break
                if hc_col is None and hc_cols:
                    hc_col = hc_cols[min(idx, len(hc_cols)-1)]

    if hdf is not None and hc_col is not None:
        h_ids = [norm_id(x) for x in hdf[sample_col].astype(str).tolist()]  # type: ignore
        h_labels = [str(x) for x in hdf[hc_col].astype(str).tolist()]  # type: ignore
        mapping = dict(zip(h_ids, h_labels))
        clabels = [mapping.get(norm_id(s), 'NA') for s in names]
    else:
        clabels = ['NA' for _ in names]

    # prepare plotting colors
    cnt = Counter(clabels)
    TOP_N = 10
    top = [c for c, _ in cnt.most_common(TOP_N) if c != 'NA']
    base_cmap = plt.get_cmap('tab20')
    colors_list = [base_cmap(i % 20) for i in range(len(top))]
    plotted_colors = [('#d0d0d0' if (v == 'NA' or v not in top) else colors_list[top.index(v)]) for v in clabels]

    if embedding_umap is not None:
        plot_emb(embedding_umap, OUT / 'umap_features_pretty.png', 'UMAP_features', names, clabels, top, colors_list, plotted_colors, hc_col, sample_counts)
        umap_data = pd.DataFrame({
            'sample_id': [str(s).lstrip('#') for s in names],
            'cluster': clabels,
            'umap1': embedding_umap[:, 0],  # type: ignore
            'umap2': embedding_umap[:, 1],  # type: ignore
        })
        umap_data.to_csv(OUT / 'umap_features_coords.csv', index=False)

    if mds_emb is not None:
        plot_emb(mds_emb, OUT / 'mds_pretty.png', 'MDS_precomputed', names, clabels, top, colors_list, plotted_colors, hc_col, sample_counts)
        mds_data = pd.DataFrame({
            'sample_id': [str(s).lstrip('#') for s in names],
            'cluster': clabels,
            'mds1': mds_emb[:, 0],  # type: ignore
            'mds2': mds_emb[:, 1],  # type: ignore
        })
        mds_data.to_csv(OUT / 'mds_coords.csv', index=False)

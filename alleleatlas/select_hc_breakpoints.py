"""Programmatic version of scripts/select_hc_breakpoints.py

Provides a function `select_and_plot_breakpoints(eval_tsv, out_prefix, top)`
that parses a hierCC eval TSV, selects candidate breakpoints, detects plateaus,
and writes silhouette + NMI images to disk.
"""
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def parse_eval_tsv(path):
    # resolve common filename variations (with/without .tsv or .gz)
    p = Path(path)
    # build candidate paths robustly: try the exact path, path + common suffixes, and replaced suffixes
    candidates = []
    candidates.append(p)
    candidates.append(Path(str(p) + '.tsv'))
    candidates.append(Path(str(p) + '.tsv.gz'))
    candidates.append(Path(str(p) + '.gz'))
    candidates.append(Path(str(p) + '.eval.tsv'))
    candidates.append(Path(str(p) + '.eval.tsv.gz'))
    # try replacing suffix with .tsv/.tsv.gz/.gz
    try:
        candidates.append(p.with_suffix('.tsv'))
    except Exception:
        pass
    try:
        candidates.append(p.with_suffix('.tsv.gz'))
    except Exception:
        pass
    try:
        candidates.append(p.with_suffix('.gz'))
    except Exception:
        pass
    # deduplicate while preserving order
    seen = set()
    unique_candidates = []
    for c in candidates:
        if str(c) not in seen:
            seen.add(str(c))
            unique_candidates.append(c)
    candidates = unique_candidates
    found = None
    for c in candidates:
        if c.exists():
            found = c
            break
    if found is None:
        raise FileNotFoundError(f"hierCC eval file not found (tried: {', '.join(str(x) for x in candidates)}). Make sure evalHCC wrote the file before calling breakpoint selection.")
    lines = found.read_text().splitlines()
    silhouette = []
    similarity_lines = []
    labels = None
    reading_nmi = False
    for ln in lines:
        if ln.startswith('#Silhouette'):
            parts = ln.strip().split('\t')
            if len(parts) >= 3:
                silhouette.append(float(parts[2]))
        if ln.startswith('#NMI'):
            parts = ln.strip().split('\t')
            if len(parts) > 1:
                labels = parts[1:]
            reading_nmi = True
            continue
        if reading_nmi:
            if ln.strip() == '':
                break
            similarity_lines.append(ln.strip().split('\t')[1:])

    silhouette = np.array(silhouette)
    if similarity_lines:
        sim = np.array([[float(x) for x in row] for row in similarity_lines])
    else:
        sim = None
    return silhouette, sim, labels


def detect_plateaus(metric, tol=1e-3, min_len=2):
    plateaus = []
    n = len(metric)
    start = 0
    i = 1
    while i < n:
        if abs(metric[i] - metric[i-1]) <= tol:
            i += 1
        else:
            if i - start >= min_len:
                plateaus.append((start, i))
            start = i
            i += 1
    if n - start >= min_len:
        plateaus.append((start, n))
    return plateaus


def select_breakpoints(silhouette, sim=None, stepwise=1,
                       sil_min=0.15, peak_prom=0.03, top_n=5):
    # use raw silhouette
    sil_smooth = silhouette.copy()
    peaks, props = find_peaks(sil_smooth, prominence=peak_prom, height=sil_min)

    n_levels = len(silhouette)
    nmi_drop = np.zeros(n_levels)
    if sim is not None:
        for i in range(n_levels-1):
            try:
                nmi_next = sim[i, i+1]
            except Exception:
                nmi_next = np.nanmean(sim[i, i+1:]) if i+1 < sim.shape[1] else 1.0
            nmi_drop[i] = 1.0 - nmi_next
        if nmi_drop.max() > 0:
            nmi_drop = nmi_drop / nmi_drop.max()
    else:
        nmi_drop = np.zeros(n_levels)

    if sil_smooth.max() > 0:
        sil_norm = (sil_smooth - sil_smooth.min()) / (sil_smooth.max() - sil_smooth.min())
    else:
        sil_norm = sil_smooth

    candidates = []
    for i in peaks:
        s_score = sil_norm[i]
        n_score = nmi_drop[i] if len(nmi_drop) > i else 0.0
        combined = 0.6 * s_score + 0.4 * n_score
        candidates.append((i, combined, s_score, n_score, sil_smooth[i]))

    candidates.sort(key=lambda x: x[1], reverse=True)
    chosen = candidates[:top_n]
    chosen_levels = [(idx*stepwise, idx, score, s, n) for idx, score, s, n, raw in chosen]
    return sil_smooth, peaks, candidates, chosen_levels


def plot_and_report(silhouette, sil_smooth, sim, peaks, candidates, chosen_levels, out_prefix, labels=None, sil_plateaus=None, nmi_plateaus=None):
    levels = np.arange(len(silhouette))
    plt.figure(figsize=(10,6))
    plt.plot(levels, silhouette, label='silhouette')

    # Shade silhouette plateaus
    sil_patch = None
    nmi_patch = None
    if sil_plateaus:
        for (a, b) in sil_plateaus:
            plt.axvspan(a, b-1, color='C2', alpha=0.12)
        sil_patch = mpatches.Patch(color='C2', alpha=0.12, label='Silhouette plateau')
    if nmi_plateaus:
        for (a, b) in nmi_plateaus:
            plt.axvspan(a, b-1, color='C3', alpha=0.08)
        nmi_patch = mpatches.Patch(color='C3', alpha=0.08, label='NMI plateau')

    # Top breakpoints
    chosen_indices = [idx for hc, idx, score, s, n in chosen_levels]
    if chosen_indices:
        marker_x = []
        marker_y = []
        marker_labels = []
        for hc, idx, score, s, n in chosen_levels:
            label = labels[idx] if (labels is not None and idx < len(labels)) else f'HC{hc}'
            plt.axvline(idx, linestyle='--', color='gray', alpha=0.6)
            marker_x.append(idx)
            marker_y.append(silhouette[idx])
            marker_labels.append(label)
        plt.scatter(marker_x, marker_y, color='red', marker='X', s=60, zorder=5, label=f'Top-{len(marker_x)} breakpoints')
        for x, y, lab in zip(marker_x, marker_y, marker_labels):
            plt.text(x, y, lab, rotation=90, va='bottom', fontsize=8)

    # show all labels on x axis (keeps previous behavior)
    n = len(silhouette)
    all_indices = list(range(n))
    all_labels = [labels[i] if (labels is not None and i < len(labels)) else str(i) for i in all_indices]
    plt.xticks(all_indices, all_labels, rotation=90, fontsize=6)
    plt.xlabel('HC level (label)')
    plt.ylabel('Silhouette')

    ax = plt.gca()
    handles, labs = ax.get_legend_handles_labels()
    extra_handles = []
    if sil_patch is not None:
        extra_handles.append(sil_patch)
    if nmi_patch is not None:
        extra_handles.append(nmi_patch)
    if extra_handles:
        plt.legend(handles=handles + extra_handles)
    else:
        plt.legend()

    plt.tight_layout()
    plt.savefig(out_prefix + '_silhouette.png', dpi=150)
    plt.close()

    if sim is not None:
        import seaborn as sns
        plt.figure(figsize=(8,6))
        sns.heatmap(sim, cmap='RdBu_r', center=0.5)
        plt.title('NMI similarity matrix')
        ax = plt.gca()
        heat_indices = list(range(sim.shape[0]))
        heat_labels = [labels[i] if (labels is not None and i < len(labels)) else str(i) for i in heat_indices]
        ax.set_xticks(heat_indices)
        ax.set_xticklabels(heat_labels, rotation=90, fontsize=6)
        ax.set_yticks(heat_indices)
        ax.set_yticklabels(heat_labels, rotation=0, fontsize=6)
        plt.tight_layout()
        plt.savefig(out_prefix + '_nmi_heatmap.png', dpi=150)
        plt.close()


def select_and_plot_breakpoints(eval_tsv, out_prefix='hiercc_breakpoints', top=5):
    silhouette, sim, labels = parse_eval_tsv(eval_tsv)
    sil_smooth, peaks, candidates, chosen_levels = select_breakpoints(silhouette, sim, top_n=top)

    sil_plateaus = detect_plateaus(sil_smooth, tol=1e-4, min_len=2)
    nmi_plateaus = []
    if sim is not None:
        n = len(silhouette)
        nmi_series = np.ones(n)
        for i in range(n-1):
            try:
                nmi_series[i] = sim[i, i+1]
            except Exception:
                nmi_series[i] = np.nanmean(sim[i, i+1:]) if i+1 < sim.shape[1] else 1.0
        nmi_series[-1] = nmi_series[-2] if n > 1 else 1.0
        nmi_diss = 1.0 - nmi_series
        nmi_plateaus = detect_plateaus(nmi_diss, tol=1e-3, min_len=2)

    plot_and_report(silhouette, sil_smooth, sim, peaks, candidates, chosen_levels, out_prefix, labels=labels, sil_plateaus=sil_plateaus, nmi_plateaus=nmi_plateaus)
    return {
        'chosen': chosen_levels,
        'sil_plateaus': sil_plateaus,
        'nmi_plateaus': nmi_plateaus,
    }


def plot_group_counts_from_hiercc(hiercc_path, eval_tsv=None, out_prefix='hiercc_group_counts', top=5):
    """Read a hierCC.HierCC.gz table and plot number of unique groups per HC level.

    If `eval_tsv` is provided, derive chosen breakpoints and plateaus from it and
    overlay them on the plot.
    Returns a dict with 'counts', 'hc_cols', and any detected breakpoints/plateaus.
    """
    import pandas as pd
    p = Path(hiercc_path)
    if not p.exists():
        raise FileNotFoundError(f'HierCC file not found: {hiercc_path}')

    # Read table (support gz compression)
    compression = 'gzip' if p.suffix == '.gz' else None
    df = pd.read_csv(p, sep='\t', compression=compression, dtype=str)

    # Heuristic: HC columns typically start with 'HC' in header. Otherwise assume all but first column
    hc_cols = [c for c in df.columns if str(c).upper().startswith('HC')]
    if not hc_cols and df.shape[1] > 1:
        hc_cols = list(df.columns[1:])

    counts = []
    for c in hc_cols:
        # treat empty strings as missing
        s = df[c].replace('', pd.NA)
        counts.append(int(s.nunique(dropna=True)))

    # If eval_tsv given, compute breakpoints and plateaus (without re-plotting the silhouette)
    chosen_levels = []
    sil_plateaus = []
    nmi_plateaus = []
    labels = None
    if eval_tsv is not None:
        try:
            silhouette, sim, labels = parse_eval_tsv(eval_tsv)
            sil_smooth, peaks, candidates, chosen_levels = select_breakpoints(silhouette, sim, top_n=top)
            sil_plateaus = detect_plateaus(sil_smooth, tol=1e-4, min_len=2)
            if sim is not None:
                n = len(silhouette)
                nmi_series = np.ones(n)
                for i in range(n-1):
                    try:
                        nmi_series[i] = sim[i, i+1]
                    except Exception:
                        nmi_series[i] = np.nanmean(sim[i, i+1:]) if i+1 < sim.shape[1] else 1.0
                nmi_series[-1] = nmi_series[-2] if n > 1 else 1.0
                nmi_plateaus = detect_plateaus(1.0 - nmi_series, tol=1e-3, min_len=2)
        except Exception:
            # don't fail plotting the counts if eval parsing fails
            chosen_levels = []
            sil_plateaus = []
            nmi_plateaus = []

    # Plot counts
    plt.figure(figsize=(10,6))
    x = np.arange(len(hc_cols))
    # reduce marker size for each data point for readability
    plt.plot(x, counts, marker='o', markersize=3, linewidth=0.8, linestyle='-', label='unique groups')
    plt.ylabel('Number of unique groups')
    plt.xlabel('HC level')
    # Reduce tick crowding: show at most `max_ticks` labels (choose evenly spaced indices)
    max_ticks = 40
    if len(hc_cols) <= max_ticks:
        tick_indices = x
        tick_labels = hc_cols
    else:
        step = max(1, len(hc_cols) // max_ticks)
        tick_indices = x[::step]
        tick_labels = [hc_cols[i] for i in tick_indices]
    plt.xticks(tick_indices, tick_labels, rotation=90, fontsize=6)

    # Overlay plateaus as shaded spans in index-space
    def _range_to_pos(a, b, labels, hc_cols):
        """Map a plateau range [a,b) in silhouette index space to positions in hc_cols.

        Returns (pos_a, pos_b) inclusive indexes for plotting axvspan(pos_a, pos_b).
        """
        # If labels available, try to map label at index a and b-1 to hc_cols positions
        if labels is not None and a < len(labels) and (b-1) < len(labels):
            la = labels[a]
            lb = labels[b-1]
            if la in hc_cols and lb in hc_cols:
                return hc_cols.index(la), hc_cols.index(lb)
        # fallback: try numeric suffix matching using helper above
        import re
        def num_from(s):
            if s is None:
                return None
            m = re.search(r"(\d+)$", str(s))
            return int(m.group(1)) if m else None

        la_num = num_from(labels[a]) if labels is not None and a < len(labels) else None
        lb_num = num_from(labels[b-1]) if labels is not None and (b-1) < len(labels) else None
        hc_nums = [num_from(c) for c in hc_cols]
        # find indices in hc_cols closest to la_num and lb_num
        def closest_index(target):
            if target is None:
                return None
            valid = [(i, n) for i, n in enumerate(hc_nums) if n is not None]
            if not valid:
                return None
            return min(valid, key=lambda pair: abs(pair[1] - target))[0]

        pa = closest_index(la_num) if la_num is not None else None
        pb = closest_index(lb_num) if lb_num is not None else None
        # If both found, return them; otherwise approximate by clamping a/b to hc_cols range
        if pa is not None and pb is not None:
            return pa, pb
        # clamp a,b into hc_cols index space
        na = int(max(0, min(len(hc_cols)-1, a)))
        nb = int(max(0, min(len(hc_cols)-1, b-1)))
        return na, nb

    if sil_plateaus:
        for (a, b) in sil_plateaus:
            pa, pb = _range_to_pos(a, b, labels, hc_cols)
            plt.axvspan(pa, pb, color='C2', alpha=0.12)
    if nmi_plateaus:
        for (a, b) in nmi_plateaus:
            pa, pb = _range_to_pos(a, b, labels, hc_cols)
            plt.axvspan(pa, pb, color='C3', alpha=0.08)

    # Helper: map silhouette index -> position in hc_cols (by matching label strings or nearest numeric match)
    def _map_idx_to_pos(idx, labels, hc_cols):
        # try direct label matching
        if labels is None:
            return idx if idx < len(hc_cols) else len(hc_cols)-1
        lbl = labels[idx] if idx < len(labels) else None
        if lbl in hc_cols:
            return hc_cols.index(lbl)
        # try numeric suffix matching
        import re
        def num_from(s):
            if s is None:
                return None
            m = re.search(r"(\d+)$", str(s))
            return int(m.group(1)) if m else None
        tgt = num_from(lbl)
        hc_nums = [num_from(c) for c in hc_cols]
        # find exact numeric match first
        if tgt is not None and tgt in hc_nums:
            return hc_nums.index(tgt)
        # otherwise find nearest numeric
        valid = [(i, n) for i, n in enumerate(hc_nums) if n is not None]
        if not valid:
            return idx if idx < len(hc_cols) else len(hc_cols)-1
        nearest = min(valid, key=lambda pair: abs(pair[1] - (tgt if tgt is not None else pair[1])))
        return nearest[0]

    # Overlay chosen breakpoints (map to hc_cols positions)
    if chosen_levels:
        marker_x = []
        marker_y = []
        for hc, idx, score, s, n in chosen_levels:
            pos = _map_idx_to_pos(idx, labels, hc_cols)
            marker_x.append(pos)
            marker_y.append(counts[pos] if pos < len(counts) else 0)
        plt.scatter(marker_x, marker_y, color='red', marker='X', s=120, label=f'Top-{len(marker_x)} breakpoints')

    # Add legend patches for plateaus
    extra_handles = []
    if sil_plateaus:
        extra_handles.append(mpatches.Patch(color='C2', alpha=0.12, label='Silhouette plateau'))
    if nmi_plateaus:
        extra_handles.append(mpatches.Patch(color='C3', alpha=0.08, label='NMI plateau'))

    ax = plt.gca()
    handles, labs = ax.get_legend_handles_labels()
    if extra_handles:
        plt.legend(handles=handles + extra_handles)
    else:
        plt.legend()

    plt.tight_layout()
    plt.savefig(out_prefix + '_group_counts.png', dpi=150)
    plt.close()

    return {
        'hc_cols': hc_cols,
        'counts': counts,
        'chosen': chosen_levels,
        'sil_plateaus': sil_plateaus,
        'nmi_plateaus': nmi_plateaus,
    }

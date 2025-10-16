import re
from pathlib import Path
from multiprocessing import Pool
import pandas as pd
from rich.console import Console

from alleleatlas.core.input import read_file_preview, detect_and_normalize_profile
from alleleatlas.cluster.pHierCC import phierCC, prepare_mat
from alleleatlas.cluster.HCCeval import evalHCC
from alleleatlas.cluster.getDistance import getDistance
from alleleatlas.select_hc_breakpoints import select_and_plot_breakpoints, plot_group_counts_from_hiercc

console = Console()

# Aliases for backward compatibility
_read_file_start = read_file_preview
adjust_input_type = detect_and_normalize_profile


def _collapse_large_profile(input_df, st_counts, outdir, target_nodes):
    """Collapse large profile to target number of clusters if needed.
    
    Returns: (normalized_path, collapse_metadata, st_counts)
    """
    from alleleatlas.cluster.collapse_profiles import auto_collapse_and_save
    
    collapse_metadata = None
    n_original_samples = sum(st_counts.values())
    
    if target_nodes is not None and input_df.shape[0] > target_nodes:
        console.print("\n[bold]Collapsing large profile[/bold]")
        console.print(f"  Original unique profiles: {input_df.shape[0]}")
        console.print(f"  Total with duplicates: {n_original_samples}")
        console.print(f"  Target clusters: {target_nodes}")
        
        temp_normalized_path = outdir + '/temp_cgmlst_profiles.normalized.tsv'
        input_df.to_csv(temp_normalized_path, sep='\t', index=False)
        
        try:
            collapsed_profile, collapse_metadata = auto_collapse_and_save(
                temp_normalized_path, outdir, target_nodes=target_nodes, nproc=4
            )
            console.print(f"[green]✓[/green] Collapsed to {collapse_metadata['n_collapsed']} clusters")
            console.print(f"  Threshold: {collapse_metadata['collapse_threshold']:.1f}")
            
            normalized_path = collapsed_profile
            st_counts = collapse_metadata['sample_counts']
            return normalized_path, collapse_metadata, st_counts
        except Exception as e:
            console.print(f"[yellow]⚠[/yellow] Collapse failed: {e}")
            console.print("[yellow]  Continuing with original profile[/yellow]")
            normalized_path = temp_normalized_path
            return normalized_path, None, st_counts
    
    # No collapsing needed
    console.print("\n[bold]Profile summary:[/bold]")
    console.print(f"  Unique profiles: {input_df.shape[0]}")
    console.print(f"  Total with duplicates: {n_original_samples}")
    
    Path(outdir).mkdir(parents=True, exist_ok=True)
    normalized_path = outdir + '/cgmlst_profiles.normalized.tsv'
    input_df.to_csv(normalized_path, sep='\t', index=False)
    return normalized_path, None, st_counts


def _compute_distances(normalized_path, nproc=4):
    """Compute both dual and pairwise distance matrices.
    
    Returns: (dist_dual, dist_p)
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


def _extract_clustering_thresholds(results, labels):
    """Extract and compute minimum clustering threshold from breakpoints and plateaus.
    
    Returns: (min_clustering_threshold, mst_drawing_threshold, bp_info)
    """
    chosen = results.get("chosen", [])
    sil_plateaus = results.get("sil_plateaus", [])
    nmi_plateaus = results.get("nmi_plateaus", [])
    
    bp_info = []
    bp_indices = []
    
    # Extract breakpoints
    for hc_level, idx, score, s, n in chosen:
        hc_label = labels[idx] if (labels is not None and idx < len(labels)) else f'HC{hc_level}'
        bp_info.append(hc_label)
        if hc_label.startswith('HC'):
            try:
                hc_value = int(hc_label[2:])
                bp_indices.append(hc_value)
            except ValueError:
                bp_indices.append(hc_level)
        else:
            bp_indices.append(hc_level)
    
    # Collect all candidate thresholds
    all_candidate_thresholds = list(bp_indices)
    
    all_plateaus = sil_plateaus + nmi_plateaus
    if all_plateaus:
        for (p_start, p_end) in all_plateaus:
            p_start_label = labels[p_start] if (labels is not None and p_start < len(labels)) else f'HC{p_start*100}'
            if p_start_label.startswith('HC'):
                try:
                    hc_value = int(p_start_label[2:])
                    all_candidate_thresholds.append(hc_value)
                except ValueError:
                    all_candidate_thresholds.append(p_start * 100)
            else:
                all_candidate_thresholds.append(p_start * 100)
    
    # Find min and second thresholds
    min_clustering_threshold = 0
    mst_drawing_threshold = None
    if all_candidate_thresholds:
        sorted_thresholds = sorted(set(all_candidate_thresholds))
        min_clustering_threshold = sorted_thresholds[0]
        if len(sorted_thresholds) > 1:
            mst_drawing_threshold = sorted_thresholds[1]
        else:
            mst_drawing_threshold = sorted_thresholds[0]
    
    return min_clustering_threshold, mst_drawing_threshold, bp_info, sil_plateaus, nmi_plateaus


def _run_breakpoint_analysis(outdir, found_eval, force):
    """Extract breakpoints and determine clustering thresholds.
    
    Parameters:
        outdir: output directory (Path)
        found_eval: path to evalHCC output file
        force: if True, recompute even if cached
        
    Returns:
        (min_clustering_threshold, mst_drawing_threshold, bp_info, sil_plateaus, nmi_plateaus)
    """
    breakpoints_plot = outdir / 'breakpoints.png'
    min_clustering_threshold = None
    mst_drawing_threshold = None
    bp_info = []
    sil_plateaus = []
    nmi_plateaus = []
    
    if found_eval is not None and (not breakpoints_plot.exists() or force):
        try:
            from alleleatlas.select_hc_breakpoints import parse_eval_tsv
            results = select_and_plot_breakpoints(found_eval, out_prefix=str(outdir / 'breakpoints'), top=5)
            silhouette, sim, labels = parse_eval_tsv(found_eval)
            
            min_clustering_threshold, mst_drawing_threshold, bp_info, sil_plateaus, nmi_plateaus = _extract_clustering_thresholds(
                results, labels
            )
            
            # Report results
            console.print(f'[bold]Suggested breakpoints:[/bold] {", ".join(bp_info)}')
            
            if sil_plateaus:
                plateau_info = []
                for (p_start, p_end) in sil_plateaus:
                    p_start_label = labels[p_start] if (labels is not None and p_start < len(labels)) else f'HC{p_start*100}'
                    p_end_label = labels[p_end-1] if (labels is not None and p_end-1 < len(labels)) else f'HC{(p_end-1)*100}'
                    plateau_info.append(f"{p_start_label}-{p_end_label}")
                console.print(f'[bold]Silhouette plateaus:[/bold] {", ".join(plateau_info)}')
            
            if nmi_plateaus:
                nmi_info = []
                for (p_start, p_end) in nmi_plateaus:
                    p_start_label = labels[p_start] if (labels is not None and p_start < len(labels)) else f'HC{p_start*100}'
                    p_end_label = labels[p_end-1] if (labels is not None and p_end-1 < len(labels)) else f'HC{(p_end-1)*100}'
                    nmi_info.append(f"{p_start_label}-{p_end_label}")
                console.print(f'[bold]NMI plateaus:[/bold] {", ".join(nmi_info)}')
            
            console.print(f'[bold]Minimum clustering threshold:[/bold] {min_clustering_threshold}')
            
        except Exception as e:
            console.print(f'[yellow]⚠[/yellow] Failed to select/plot breakpoints: {e}')
    elif breakpoints_plot.exists() and not force:
        console.print(f'[yellow]↷[/yellow] Using cached breakpoints: {breakpoints_plot}')
    else:
        console.print('[yellow]⚠[/yellow] hierCC eval file not found; skipping breakpoint selection.')
    
    return min_clustering_threshold, mst_drawing_threshold, bp_info, sil_plateaus, nmi_plateaus


def _run_group_counts(outdir, hiercc_path, found_eval, force):
    """Generate group-counts visualization from HierCC clustering.
    
    Parameters:
        outdir: output directory (Path)
        hiercc_path: path to hierCC.HierCC.gz file
        found_eval: path to evalHCC output file
        force: if True, recompute even if cached
    """
    gc_plot = outdir / 'group_counts_summary.png'
    if hiercc_path.exists() and found_eval is not None and (not gc_plot.exists() or force):
        try:
            gc_results = plot_group_counts_from_hiercc(str(hiercc_path), eval_tsv=found_eval, out_prefix=str(outdir / 'group_counts'), top=5)
            gc_chosen = gc_results.get("chosen", [])
            
            if gc_chosen:
                from alleleatlas.select_hc_breakpoints import parse_eval_tsv
                silhouette, sim, labels = parse_eval_tsv(found_eval)
                gc_bp_info = []
                for hc_level, idx, score, s, n in gc_chosen:
                    hc_label = labels[idx] if (labels is not None and idx < len(labels)) else f'HC{hc_level}'
                    gc_bp_info.append(hc_label)
                console.print(f'[green]✓[/green] Group-counts breakpoints: {", ".join(gc_bp_info)}')
            else:
                console.print('[green]✓[/green] Group-counts plot written')
        except Exception as e:
            console.print(f'[yellow]⚠[/yellow] Failed to produce group-counts plot: {e}')
    elif gc_plot.exists() and not force:
        console.print(f'[yellow]↷[/yellow] Using cached group-counts plot: {gc_plot}')


def _run_mst_visualization(outdir, normalized_path, min_threshold, mst_threshold, collapse_metadata, dist_dual, force):
    """Generate MST network visualization with group coloring.
    
    Parameters:
        outdir: output directory (Path)
        normalized_path: path to normalized profile
        min_threshold: minimum clustering threshold for node collapse
        mst_threshold: line threshold for edges
        collapse_metadata: metadata about collapsed samples
        dist_dual: precomputed distance matrix
        force: if True, recompute even if cached
    """
    mst_plot = outdir / 'mst_visualization.png'
    if mst_threshold is not None and (not mst_plot.exists() or force):
        # NOTE: MST drawing functionality not yet implemented
        # Placeholder for future Phase 3 implementation (visualization/mst.py)
        console.print('[yellow]⚠[/yellow] MST drawing not yet implemented (Phase 3 refactoring)')
    elif mst_plot.exists() and not force:
        console.print('[yellow]↷[/yellow] Using cached MST visualization')
    else:
        console.print('[yellow]⚠[/yellow] No clustering thresholds determined; skipping MST.')


def _run_umap_embeddings(outdir, normalized_path, force):
    """Generate UMAP dimensionality reduction embeddings.
    
    Parameters:
        outdir: output directory (Path)
        normalized_path: path to normalized profile
        force: if True, recompute even if cached
    """
    umap_plot = outdir / 'umap_embeddings.png'
    if not umap_plot.exists() or force:
        try:
            from alleleatlas.umap_from_distance import run as umap_run
            console.print('[bold]Computing UMAP embeddings...[/bold]')
            umap_run(str(normalized_path), str(outdir), nproc=2)
            console.print('[green]✓[/green] UMAP embeddings complete')
        except Exception as e:
            console.print(f'[yellow]⚠[/yellow] Failed to compute UMAP: {e}')
    else:
        console.print('[yellow]↷[/yellow] Using cached UMAP embeddings')


def main(cgmlst_profiles, outdir, target_nodes=None, force=False):
    """Main pipeline orchestrator for cgMLST clustering analysis.
    
    Eight-step pipeline:
    1. Load and normalize input (auto-detect format)
    2. Collapse samples if needed (for large datasets)
    3. Compute pairwise distance matrices (cached)
    4. Hierarchical clustering (HierCC)
    5. Stability evaluation (breakpoints, plateaus)
    6. Group-counts visualization
    7. UMAP embeddings
    8. MST network visualization
    
    Parameters:
        cgmlst_profiles (str): path to cgMLST profile file
        outdir (str): output directory
        target_nodes (int): target number of clusters for collapsing (auto-detect if None)
        force (bool): if True, recompute all steps even if cached (default: False)
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # ========================================================================
    # STEP 1: Load and normalize input
    # ========================================================================
    input_df, st_counts = adjust_input_type(cgmlst_profiles)
    if target_nodes is None:
        target_nodes = 10000 if input_df.shape[0] > 10000 else None
    
    # ========================================================================
    # STEP 2: Collapse if needed (or use cached if exists and not forced)
    # ========================================================================
    normalized_path = outdir / 'cgmlst_profiles.normalized.tsv'
    if normalized_path.exists() and not force:
        console.print(f'[yellow]↷[/yellow] Using cached normalized profile: {normalized_path}')
    else:
        normalized_path, collapse_metadata, st_counts = _collapse_large_profile(
            input_df, st_counts, str(outdir), target_nodes
        )
    console.print(f'Normalized cgMLST profile: {normalized_path}')
    
    # ========================================================================
    # STEP 3: Compute distances (or use cached if exists and not forced)
    # ========================================================================
    dist_dual, dist_p = _compute_distances(str(normalized_path), nproc=4)
    
    # ========================================================================
    # STEP 4: Run clustering and evaluation (or skip if cached and not forced)
    # ========================================================================
    hiercc_output = outdir / 'hierCC.HierCC.gz'
    if hiercc_output.exists() and not force:
        console.print(f'[yellow]↷[/yellow] Using cached clustering: {hiercc_output}')
    else:
        console.print('[bold]Running hierarchical clustering and evaluation...[/bold]')
        phierCC(str(normalized_path), str(outdir / 'hierCC'), '', 4, 0.03, precomputed_dist=dist_dual)
        evalHCC(str(normalized_path), str(outdir / 'hierCC.HierCC.gz'), str(outdir / 'hierCC.eval'), 100, 4, precomputed_dist=dist_p)
    
    # ========================================================================
    # STEP 5: Extract breakpoints and determine thresholds (or skip if cached)
    # ========================================================================
    candidates = [str(outdir / f) for f in ['hierCC.eval', 'hierCC.eval.tsv', 'hierCC.eval.tsv.gz', 'hierCC.eval.gz']]
    found_eval = None
    for c in candidates:
        if Path(c).exists():
            found_eval = c
            break
    
    collapse_metadata = None
    min_clustering_threshold, mst_drawing_threshold, bp_info, sil_plateaus, nmi_plateaus = _run_breakpoint_analysis(
        outdir, found_eval, force
    )
    
    # ========================================================================
    # STEP 6: Generate group-counts plot (or skip if cached)
    # ========================================================================
    _run_group_counts(outdir, outdir / 'hierCC.HierCC.gz', found_eval, force)
    
    # ========================================================================
    # STEP 7: Draw MST visualization (or skip if cached)
    # ========================================================================
    _run_mst_visualization(outdir, normalized_path, min_clustering_threshold, mst_drawing_threshold, collapse_metadata, dist_dual, force)
    
    # ========================================================================
    # STEP 8: Generate UMAP embeddings (or skip if cached)
    # ========================================================================
    _run_umap_embeddings(outdir, normalized_path, force)
    

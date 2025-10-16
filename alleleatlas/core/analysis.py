"""Analysis functions for hierarchical clustering results.

Includes threshold extraction and plateau detection from evaluation metrics.
"""


def extract_clustering_thresholds(results, labels):
    """Extract and compute minimum clustering threshold from breakpoints and plateaus.
    
    Parameters:
        results (dict): Breakpoint analysis results containing:
            - "chosen": list of chosen breakpoints
            - "sil_plateaus": silhouette plateaus as (start, end) tuples
            - "nmi_plateaus": NMI plateaus as (start, end) tuples
        labels (list): HC threshold labels corresponding to indices
        
    Returns:
        tuple: (min_clustering_threshold, mst_drawing_threshold, bp_info, sil_plateaus, nmi_plateaus)
            - min_clustering_threshold: Recommended minimum threshold for clustering
            - mst_drawing_threshold: Secondary threshold for MST visualization
            - bp_info: List of breakpoint labels for reporting
            - sil_plateaus: Silhouette plateaus
            - nmi_plateaus: NMI plateaus
            
    Notes:
        - Extracts HC thresholds from chosen breakpoints
        - Identifies minimum candidate threshold from all sources
        - Selects second-lowest threshold for MST visualization if available
    """
    chosen = results.get("chosen", [])
    sil_plateaus = results.get("sil_plateaus", [])
    nmi_plateaus = results.get("nmi_plateaus", [])
    
    bp_info = []
    bp_indices = []
    
    # Extract breakpoints
    for hc_level, idx, *_ in chosen:
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
        for (p_start, _) in all_plateaus:
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


# Backward compatibility alias
_extract_clustering_thresholds = extract_clustering_thresholds


__all__ = ['extract_clustering_thresholds']

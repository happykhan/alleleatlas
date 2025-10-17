"""Deduplication and pre-collapse at exact match (HC0) level.

For large datasets (>10k profiles), deduplicates exact matches before
distance computation to achieve 10x+ speedup on datasets with many duplicates.
"""

from rich.console import Console

console = Console()


def deduplicate_at_exact_match(input_df, st_counts):
    """Pre-collapse profiles at exact match (HC0) level.
    
    For large datasets with many exact duplicates, this dramatically reduces
    the distance computation burden from O(n²) → O(m²) where m << n.
    
    Parameters:
        input_df (pd.DataFrame): Normalized cgMLST profiles
        st_counts (dict): Sample counts per profile
        
    Returns:
        tuple: (deduplicated_df, dedup_mapping, st_counts_updated)
            - deduplicated_df: Profiles with exact duplicates removed
            - dedup_mapping: Dict mapping deduplicated indices to original profiles
            - st_counts_updated: Updated sample counts
            
    Notes:
        - Exact match = all loci identical (HC0)
        - Keeps first occurrence of each unique profile
        - Updates st_counts to sum duplicates
        - Fast: O(n) instead of distance computation O(n²)
    """
    n_original = len(input_df)
    
    # Get profile columns (all except ID/ST columns)
    profile_cols = input_df.columns[1:].tolist()
    
    # Find duplicate profiles (exact matches)
    # Use .drop_duplicates() to keep only first occurrence
    dedup_df = input_df.drop_duplicates(subset=profile_cols, keep='first').reset_index(drop=True)
    
    n_deduplicated = len(dedup_df)
    n_removed = n_original - n_deduplicated
    
    if n_removed > 0:
        reduction_pct = 100 * n_removed / n_original
        console.print('[bold]Pre-collapse at exact match (HC0)[/bold]')
        console.print(f'  Original profiles: {n_original:,}')
        console.print(f'  Deduplicated: {n_deduplicated:,}')
        console.print(f'  Removed: {n_removed:,} ({reduction_pct:.1f}%)')
        
        # Calculate speedup
        speedup_factor = (n_original / n_deduplicated) ** 2
        console.print(f'  [green]Expected distance computation speedup: {speedup_factor:.1f}x[/green]')
        
        # Update st_counts: sum samples for duplicate profiles
        # For profiles that were removed, transfer their samples to the kept one
        dedup_st_counts = {}
        
        for orig_idx in range(n_original):
            orig_profile = input_df.iloc[orig_idx]
            
            # Find if this profile exists in deduplicated set
            matches = (dedup_df[profile_cols] == orig_profile[profile_cols]).all(axis=1)
            
            if matches.any():
                kept_idx = dedup_df[matches].index[0]
                
                # Sum counts for this profile
                if kept_idx not in dedup_st_counts:
                    dedup_st_counts[kept_idx] = 0
                    
                # Add sample count from original
                st_id = orig_profile.iloc[0]  # Assuming first column is ST/ID
                if st_id in st_counts:
                    dedup_st_counts[kept_idx] += st_counts[st_id]
                else:
                    dedup_st_counts[kept_idx] += 1
        
        # Rebuild st_counts dict with deduplicated IDs
        dedup_st_counts_final = {}
        for new_idx, row in dedup_df.iterrows():
            st_id = row.iloc[0]
            dedup_st_counts_final[st_id] = dedup_st_counts.get(new_idx, 1)
        
        return dedup_df, n_removed, dedup_st_counts_final
    else:
        console.print('[yellow]No exact duplicates found (HC0 level)[/yellow]')
        return input_df, 0, st_counts


def should_deduplicate(n_profiles, force_dedup=False):
    """Determine if pre-collapse deduplication should be applied.
    
    Parameters:
        n_profiles (int): Number of profiles in dataset
        force_dedup (bool): Force deduplication regardless of size
        
    Returns:
        bool: True if deduplication should be applied
    """
    # Apply for large datasets (>10k profiles) or if forced
    return n_profiles > 10000 or force_dedup


__all__ = [
    'deduplicate_at_exact_match',
    'should_deduplicate',
]

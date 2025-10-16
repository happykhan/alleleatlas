"""Shared visualization utilities for UMAP and other plotting functions."""

import re


def num_from(lbl):
    """Extract the numeric suffix from a label.
    
    Args:
        lbl: Label string or None
        
    Returns:
        Integer numeric suffix if found, None otherwise
        
    Examples:
        >>> num_from("cluster_5")
        5
        >>> num_from("ST123")
        123
        >>> num_from(None)
        None
    """
    if lbl is None:
        return None
    m = re.search(r"(\d+)$", str(lbl))
    return int(m.group(1)) if m else None


def norm_id(x):
    """Normalize an identifier by stripping whitespace and leading hashes.
    
    Args:
        x: Identifier to normalize
        
    Returns:
        Normalized string, empty string if input is None
        
    Examples:
        >>> norm_id("  #sample123  ")
        'sample123'
        >>> norm_id(None)
        ''
    """
    if x is None:
        return ''
    return str(x).strip().lstrip('#')


__all__ = ['num_from', 'norm_id']

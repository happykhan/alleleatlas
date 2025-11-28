"""Compatibility shim exposing `approx` modules by re-exporting
the implementations present under `alleleatlas`.
"""

__all__ = [
    'convert_to_parquet',
    'build_knn',
    'analyze_knn',
    'mst_and_plot',
    'backends',
    'analysis',
    'validation',
]

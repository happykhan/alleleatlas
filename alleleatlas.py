#!/usr/bin/env python3
"""AlleleAtlas: Clustering and visualization of cgMLST profiles.

Simple CLI entry point that orchestrates the full clustering pipeline.
"""

import argparse
import sys
from pathlib import Path

from rich.console import Console
from alleleatlas.main import run_pipeline, ClusteringConfig

console = Console()


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="AlleleAtlas: cgMLST clustering and visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python alleleatlas.py cgmlst_data/profiles.csv output_dir/
  python alleleatlas.py data.csv output/ --eps 1.5 --minpts 5
  python alleleatlas.py data.csv output/ --force
        """,
    )

    parser.add_argument(
        "input",
        type=str,
        help="Input cgMLST profile file (CSV/TSV with ST column + gene columns)",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Output directory for results",
    )

    parser.add_argument(
        "--eps",
        type=float,
        default=None,
        help="DBSCAN eps parameter (default: auto-computed)",
    )
    parser.add_argument(
        "--minpts",
        type=int,
        default=3,
        help="DBSCAN minPts parameter (default: 3)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force recomputation of distance matrices (ignore cache)",
    )
    parser.add_argument(
        "--nproc",
        type=int,
        default=2,
        help="Number of processes for parallel computation (default: 2)",
    )

    args = parser.parse_args()

    # Create config
    config = ClusteringConfig(
        cgmlst_profiles=args.input,
        outdir=args.output,
        dbscan_eps=args.eps,
        dbscan_min_samples=args.minpts,
        force_recompute=args.force,
        nproc=args.nproc,
    )

    # Run pipeline
    try:
        results = run_pipeline(config)
        console.print("\n✓ Pipeline finished successfully!")
        return 0
    except Exception as e:
        console.print(f"\n✗ Pipeline failed: {e}", style="bold red")
        return 1


if __name__ == "__main__":
    sys.exit(main())
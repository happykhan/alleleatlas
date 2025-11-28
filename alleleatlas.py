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
        "--k",
        type=int,
        default=50,
        help="Number of neighbors for k-NN (default: 50)",
    )
    parser.add_argument(
        "--backend",
        choices=["exact", "usearch"],
        default="exact",
        help="Distance backend to use",
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

    # Create config for the new approx-based pipeline (thin wrapper)
    config = ClusteringConfig(
        cgmlst_profiles=args.input,
        outdir=args.output,
        k=args.k,
        backend=args.backend,
        force_recompute=args.force,
        nproc=args.nproc,
    )

    # Run pipeline
    try:
        results = run_pipeline(config)
        console.print("\n✓ Pipeline finished successfully!")
        return 0
    except Exception as e:
        # Escape markup characters in exception message to avoid Rich parsing errors
        error_msg = str(e).replace('[', '\\[').replace(']', '\\]')
        console.print(f"\n✗ Pipeline failed: {error_msg}", style="bold red")
        return 1


if __name__ == "__main__":
    sys.exit(main())

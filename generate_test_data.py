#!/usr/bin/env python3
"""Generate synthetic cgMLST test data with known cluster groups.

This script creates small, simple datasets with preknown cluster structures
for testing and validation of the AlleleAtlas pipeline.
"""

import numpy as np
import pandas as pd
from pathlib import Path
import argparse
from typing import Tuple, List


def generate_cluster_profiles(
    cluster_id: int,
    cluster_size: int,
    n_genes: int,
    n_alleles: int,
    base_profile=None,
    mutation_rate: float = 0.05,
    seed=None,
) -> Tuple[np.ndarray, List[str]]:
    """Generate profiles for a single cluster with controlled mutation rate.
    
    Parameters:
        cluster_id: Cluster number (for naming)
        cluster_size: Number of profiles in this cluster
        n_genes: Number of genes (loci)
        n_alleles: Number of possible alleles per gene
        base_profile: Starting profile (if None, random)
        mutation_rate: Fraction of genes to mutate per profile
        seed: Random seed for reproducibility
        
    Returns:
        profiles: Array of shape (cluster_size, n_genes)
        sample_names: List of sample IDs
    """
    if seed is not None:
        np.random.seed(seed + cluster_id)
    
    # Generate base profile if not provided
    if base_profile is None:
        base_profile = np.random.randint(1, n_alleles + 1, n_genes)
    else:
        base_profile = base_profile.copy()
    
    profiles = []
    sample_names = []
    
    for i in range(cluster_size):
        # Start from base profile
        profile = base_profile.copy()
        
        # Randomly mutate some genes
        n_mutations = max(1, int(n_genes * mutation_rate))
        mutation_positions = np.random.choice(n_genes, n_mutations, replace=False)
        
        for pos in mutation_positions:
            # Change to different allele
            new_allele = np.random.randint(1, n_alleles + 1)
            while new_allele == profile[pos]:  # Ensure it's different
                new_allele = np.random.randint(1, n_alleles + 1)
            profile[pos] = new_allele
        
        profiles.append(profile)
        sample_names.append(f"C{cluster_id}_S{i:03d}")  # e.g., C0_S001
    
    return np.array(profiles), sample_names


def generate_synthetic_data(
    n_clusters: int = 3,
    cluster_size: int = 20,
    n_genes: int = 50,
    n_alleles: int = 100,
    mutation_rate: float = 0.05,
    seed: int = 42,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    """Generate synthetic cgMLST data with known clusters.
    
    Parameters:
        n_clusters: Number of clusters
        cluster_size: Samples per cluster
        n_genes: Number of genes/loci
        n_alleles: Number of possible alleles per gene
        mutation_rate: Fraction of genes mutated per profile relative to cluster base
        seed: Random seed
        
    Returns:
        df: DataFrame with columns: ST (sequence type), gene_1, gene_2, ..., gene_N
        cluster_labels: True cluster assignment for each sample
        sample_ids: Sample names
    """
    np.random.seed(seed)
    
    all_profiles = []
    all_sample_names = []
    cluster_labels = []
    
    for cluster_id in range(n_clusters):
        # Generate base profile for this cluster
        base_profile = np.random.randint(1, n_alleles + 1, n_genes)
        
        # Generate profiles in this cluster
        profiles, sample_names = generate_cluster_profiles(
            cluster_id=cluster_id,
            cluster_size=cluster_size,
            n_genes=n_genes,
            n_alleles=n_alleles,
            base_profile=base_profile,
            mutation_rate=mutation_rate,
            seed=seed,
        )
        
        all_profiles.append(profiles)
        all_sample_names.extend(sample_names)
        cluster_labels.extend([cluster_id] * len(sample_names))
    
    # Combine all profiles
    all_profiles = np.vstack(all_profiles)
    
    # Create DataFrame
    gene_names = [f"gene_{i:02d}" for i in range(n_genes)]
    df = pd.DataFrame(all_profiles, columns=gene_names)
    df.insert(0, 'ST', all_sample_names)
    
    return df, cluster_labels, all_sample_names


def generate_challenging_data(
    seed: int = 42,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    """Generate challenging dataset: clusters with varying separation.
    
    Creates:
    - Cluster A: Tight, well-separated
    - Cluster B: Medium separation
    - Cluster C: Loose, overlapping with B
    
    Returns:
        df: DataFrame
        cluster_labels: True labels
        sample_ids: Sample names
    """
    np.random.seed(seed)
    
    n_genes = 100
    n_alleles = 200
    
    all_profiles = []
    all_sample_names = []
    cluster_labels = []
    
    # Cluster A: Tight (low mutation rate)
    base_a = np.random.randint(1, n_alleles + 1, n_genes)
    profiles_a, names_a = generate_cluster_profiles(
        0, 15, n_genes, n_alleles, base_a, mutation_rate=0.02, seed=seed
    )
    all_profiles.append(profiles_a)
    all_sample_names.extend(names_a)
    cluster_labels.extend([0] * len(names_a))
    
    # Cluster B: Medium (medium mutation rate)
    base_b = base_a.copy()
    base_b[np.random.choice(n_genes, 20, replace=False)] = np.random.randint(
        1, n_alleles + 1, 20
    )
    profiles_b, names_b = generate_cluster_profiles(
        1, 15, n_genes, n_alleles, base_b, mutation_rate=0.05, seed=seed
    )
    all_profiles.append(profiles_b)
    all_sample_names.extend(names_b)
    cluster_labels.extend([1] * len(names_b))
    
    # Cluster C: Loose (high mutation rate, overlaps with B)
    base_c = base_b.copy()
    base_c[np.random.choice(n_genes, 10, replace=False)] = np.random.randint(
        1, n_alleles + 1, 10
    )
    profiles_c, names_c = generate_cluster_profiles(
        2, 15, n_genes, n_alleles, base_c, mutation_rate=0.10, seed=seed
    )
    all_profiles.append(profiles_c)
    all_sample_names.extend(names_c)
    cluster_labels.extend([2] * len(names_c))
    
    all_profiles = np.vstack(all_profiles)
    gene_names = [f"gene_{i:03d}" for i in range(n_genes)]
    df = pd.DataFrame(all_profiles, columns=gene_names)
    df.insert(0, 'ST', all_sample_names)
    
    return df, cluster_labels, all_sample_names


def save_data(
    df: pd.DataFrame,
    output_path: Path,
    compress: bool = True,
) -> None:
    """Save dataset in cgMLST format.
    
    Parameters:
        df: DataFrame with profiles
        output_path: Path to save
        compress: Save as .gz if True
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if compress:
        # Save as gzipped CSV
        df.to_csv(output_path, sep='\t', index=False)
        print(f"✓ Saved to {output_path}")
    else:
        df.to_csv(output_path, sep='\t', index=False)
        print(f"✓ Saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic cgMLST test data with known cluster groups"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="cgmlst_test_synthetic.csv",
        help="Output filename (default: cgmlst_test_synthetic.csv)"
    )
    parser.add_argument(
        "-n", "--n-clusters",
        type=int,
        default=3,
        help="Number of clusters (default: 3)"
    )
    parser.add_argument(
        "-s", "--cluster-size",
        type=int,
        default=20,
        help="Samples per cluster (default: 20)"
    )
    parser.add_argument(
        "-g", "--n-genes",
        type=int,
        default=50,
        help="Number of genes/loci (default: 50)"
    )
    parser.add_argument(
        "-a", "--n-alleles",
        type=int,
        default=100,
        help="Number of possible alleles per gene (default: 100)"
    )
    parser.add_argument(
        "-m", "--mutation-rate",
        type=float,
        default=0.05,
        help="Mutation rate per cluster (default: 0.05 = 5%%)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)"
    )
    parser.add_argument(
        "--challenging",
        action="store_true",
        help="Generate challenging dataset (varying cluster separation)"
    )
    parser.add_argument(
        "--no-compress",
        action="store_true",
        help="Don't compress output"
    )
    parser.add_argument(
        "--meta",
        action="store_true",
        help="Also save metadata with true cluster assignments"
    )
    
    args = parser.parse_args()
    
    print("\n" + "="*60)
    print("Synthetic cgMLST Test Data Generator")
    print("="*60)
    
    # Generate data
    if args.challenging:
        print("\n→ Generating CHALLENGING dataset...")
        print("  Cluster A: Tight (well-separated)")
        print("  Cluster B: Medium (moderate overlap)")
        print("  Cluster C: Loose (overlapping)")
        df, true_labels, sample_ids = generate_challenging_data(seed=args.seed)
    else:
        print(f"\n→ Generating {args.n_clusters} clusters")
        print(f"  Cluster size: {args.cluster_size} samples each")
        print(f"  Genes: {args.n_genes}")
        print(f"  Alleles per gene: {args.n_alleles}")
        print(f"  Mutation rate: {args.mutation_rate*100:.1f}%")
        df, true_labels, sample_ids = generate_synthetic_data(
            n_clusters=args.n_clusters,
            cluster_size=args.cluster_size,
            n_genes=args.n_genes,
            n_alleles=args.n_alleles,
            mutation_rate=args.mutation_rate,
            seed=args.seed,
        )
    
    # Save main data
    output_path = Path(args.output)
    save_data(df, output_path, compress=not args.no_compress)
    
    # Print summary
    print("\n[Dataset Summary]")
    print(f"  Samples: {len(df)}")
    print(f"  Genes: {len(df.columns) - 1}")
    print(f"  Unique clusters: {len(set(true_labels))}")
    
    # Count per cluster
    for cluster_id in sorted(set(true_labels)):
        count = true_labels.count(cluster_id)
        print(f"    Cluster {cluster_id}: {count} samples")
    
    # Save metadata if requested
    if args.meta:
        meta_path = output_path.parent / f"{output_path.stem}.metadata.txt"
        
        meta_df = pd.DataFrame({
            'Sample': sample_ids,
            'True_Cluster': true_labels,
        })
        meta_df.to_csv(meta_path, sep='\t', index=False)
        print(f"\n✓ Metadata saved to {meta_path}")
        print("\nMetadata content:")
        print(meta_df.to_string(index=False))
    
    print(f"\n→ Ready to run: python alleleatlas.py {output_path} output_dir")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()

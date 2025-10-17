"""AlleleAtlas Refactored Pipeline - Main Orchestrator (v2)

Nine-step clustering and visualization pipeline:
1. Input & Cleanup
2. Local Clustering (DBSCAN)
3. Global Clustering (UPGMA) + Dendrogram
4. Thresholding & Analysis
5. Statistical Validation (Permutation Tests)
6. Collapsed Node Analysis  
7. Network Visualization
8. Output Summary
9. Clustering Visualizations (MDS, Silhouette, Metrics)
"""
import os 
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any
import numpy as np
from rich.console import Console
import pandas as pd

from alleleatlas.core.input import detect_and_normalize_profile
from alleleatlas.core.distances import compute_distance_matrices
from alleleatlas.cluster.dbscan_local import run_dbscan_clustering
from alleleatlas.cluster.upgma_global import run_upgma_clustering
from alleleatlas.cluster.permutation_tests import run_permutation_tests
from alleleatlas.visualization.network_graph import run_network_visualization
from alleleatlas.visualization.clustering_metrics import (
    plot_mds_with_clusters,
    plot_silhouette_scores,
    plot_clustering_metrics,
    plot_hcceval_style,
    plot_hiercc_group_sizes,
)

from alleleatlas.cluster.myHierCC import my_phierCC
from alleleatlas.cluster.myHierCC import evalHCC

console = Console()


@dataclass
class ClusteringConfig:
    """Configuration for the clustering pipeline."""
    
    # Required
    cgmlst_profiles: str
    outdir: str
    
    # DBSCAN
    dbscan_eps: Optional[float] = None  # If None, auto-detect
    dbscan_min_samples: int = 5
    
    # UPGMA
    upgma_method: str = 'single'
    
    # Thresholding
    min_cluster_size: int = 2
    
    # Permutation testing
    n_permutations: int = 100
    n_perm_procs: int = 4
    run_permutation_tests_flag: bool = True
    
    # Network visualization
    layout_method: str = 'spring'  # 'spring', 'spectral', 'equidistant', 'circular'
    edge_threshold: Optional[float] = None  # If None, use median
    edge_percentile: Optional[float] = None  # E.g., 75 for 75th percentile
    
    # General
    force_recompute: bool = False
    nproc: int = 4
    verbose: bool = True


class ClusteringPipeline:
    """Main orchestrator for the 8-step clustering pipeline."""
    
    def __init__(self, config: ClusteringConfig):
        self.config = config
        self.outdir = Path(config.outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        
        # Result storage - typed as Any to avoid Pylance inference issues
        self.results: Dict[str, Any] = {
            'step_1_input': {},
            'step_2_dbscan': {},
            'step_3_upgma': {},
            'step_4_thresholding': {},
            'step_5_permutation': {},
            'step_6_nodes': {},
            'step_7_network': {},
            'step_8_summary': {},
        }
    
    def run(self):
        """Execute the full 8-step pipeline."""
        console.print("[bold cyan]═══════════════════════════════════════════════════════════[/bold cyan]")
        console.print("[bold cyan]AlleleAtlas Refactored Pipeline - 8-Step Clustering[/bold cyan]")
        console.print("[bold cyan]═══════════════════════════════════════════════════════════[/bold cyan]\n")
        
        try:
            self._step_1_input_cleanup()
            self.hiercc()
            self._step_2_dbscan_local_clustering()
            self._step_3_upgma_global_clustering()
            self._step_4_thresholding_analysis()
            if self.config.run_permutation_tests_flag:
                self._step_5_permutation_testing()
            self._step_7_network_visualization()
            self._step_8_summary()
            self._step_9_clustering_visualizations()
            
            console.print("\n[bold green]✓[/bold green] Pipeline completed successfully!")
            console.print(f"[bold]Output directory:[/bold] {self.outdir}")
            
            return self.results
        
        except Exception as e:
            console.print(f"\n[bold red]✗[/bold red] Pipeline failed: {e}")
            import traceback
            traceback.print_exc()
            raise
    
    def _step_1_input_cleanup(self):
        """Step 1: Load and normalize input."""
        console.print("\n[bold]STEP 1: Input & Cleanup[/bold]")
        console.print("─" * 60)
        
        # Load profile
        input_df, st_counts = detect_and_normalize_profile(self.config.cgmlst_profiles)
        console.print(f"  Loaded {len(input_df)} unique profiles")
        console.print(f"  Total samples (with duplicates): {sum(st_counts.values())}")
        
        # Save normalized profile
        # TODO: This should be set to be optional output (it's only for debugging)
        normalized_path = self.outdir / 'cgmlst_profiles.normalized.tsv.gz'
        if not normalized_path.exists() or self.config.force_recompute:
            input_df.to_csv(str(normalized_path), sep='\t', index=False, compression='gzip')
            console.print(f"  Saved normalized profile to {normalized_path}")
        # Check for cached distance matrices
        dist_dual_cache = self.outdir / 'distance_matrix.dual.npy'
        dist_p_cache = self.outdir / 'distance_matrix.pairwise.npy'
        
        # Compute distance matrices
        console.print("\n  Computing distance matrices...")
        dist_dual, dist_p, mat, names = compute_distance_matrices(str(normalized_path), nproc=self.config.nproc, allowed_missing=0.05)
        
        # Check available disk space before caching
        total_cache_size = dist_dual.nbytes + dist_p.nbytes + mat.nbytes + names.nbytes
        available_space = shutil.disk_usage(str(self.outdir)).free
        
        if available_space > total_cache_size * 1.5:  # Need 50% extra margin
            # Save to cache
            console.print("  Saving distance matrices to cache...")
            try:
                np.save(str(dist_dual_cache), dist_dual)
                np.save(str(dist_p_cache), dist_p)
                np.save(self.outdir / 'matrix.npy', mat)
                np.save(self.outdir / 'names.npy', names)
                console.print(f"  ✓ Cached {total_cache_size / 1e9:.1f} GB")
            except OSError as e:
                console.print(f"  ⚠ Could not cache matrices (disk full): {e}")
                # Continue without cache - matrices already in memory
        else:
            console.print(f"  ⚠ Insufficient disk space for cache ({available_space / 1e9:.1f} GB available, {total_cache_size / 1e9:.1f} GB needed)")
            console.print("  → Skipping cache, using in-memory matrices only")
        
        # Save human-readable compressed versions
        try:            
            # Save pairwise as compressed CSV (handle both 2D and 3D)
            dist_p_csv = self.outdir / 'distance_matrix.pairwise.csv.gz'
            if dist_p.ndim == 3:
                dist_p_rows = [row.flatten().tolist() for row in dist_p]
                pd.DataFrame(dist_p_rows).to_csv(dist_p_csv, compression='gzip', index=False)
            else:
                pd.DataFrame(dist_p).to_csv(dist_p_csv, compression='gzip', index=False)
            
            # Save dual (lower triangle) - convert to list of lists for CSV
            dist_dual_csv = self.outdir / 'distance_matrix.dual.csv.gz'
            if dist_dual.ndim == 3:
                dist_dual_rows = [row.flatten().tolist() for row in dist_dual]
            else:
                dist_dual_rows = [row.flatten().tolist() for row in dist_dual.reshape(dist_dual.shape[0], -1)]
            pd.DataFrame(dist_dual_rows).to_csv(dist_dual_csv, compression='gzip', index=False)
            console.print("  ✓ Saved compressed readable versions (.csv.gz)")
        except Exception as e:  # noqa: BLE001
            console.print(f"  ⚠ Could not save readable versions: {e}")
        
        self.results['step_1_input'] = {
            'input_df': input_df,
            'st_counts': st_counts,
            'distance_matrix': dist_p,
            'matrix_raw': mat,
            'matrix_names': names,
            'distance_dual': dist_dual,
            'distance_dual_path': str(dist_dual_cache),
            'distance_pairwise_path': str(dist_p_cache),
            'normalized_path': str(normalized_path),
        }
    


    def hiercc(self):
        """Run HierCC clustering and HierCC evaluation."""#
        console.print("\n[bold]STEP 1.5: HierCC Clustering & Evaluation[/bold]")
        console.print("─" * 60)
        dist_dual_matrix = self.results['step_1_input']['distance_dual']
        mat = self.results['step_1_input']['matrix_raw']
        names = self.results['step_1_input']['matrix_names']
        n_proc = self.config.nproc
        profiles = self.config.cgmlst_profiles
        console.print(f"  Running HierCC on profiles: {profiles}")
        output = os.path.join(self.outdir, 'hiercc')
        res = my_phierCC(mat, names, dist_dual_matrix, output, n_proc)
        # Evaluate HierCC results
        output_prefix = output + '_HCCeval'
        stepwise = 10 
        if len(mat) < 1000:
            stepwise = 1
        silhouette, similarity = evalHCC(mat, names, res, output_prefix, stepwise=stepwise, n_proc=n_proc)
        
        # Plot group sizes at each level
        group_count_plot, group_counts = plot_hiercc_group_sizes(res, stepwise, str(self.outdir))
        
        self.results['step_1.5_hiercc'] = {
            'silhouette': silhouette,
            'similarity': similarity,
            'hiercc_groups': res,
            'hiercc_output_prefix': output_prefix,
            'stepwise': stepwise,
            'group_count_plot': group_count_plot,
            'group_counts': group_counts
        }

    def _step_2_dbscan_local_clustering(self):
        """Step 2: Local clustering via DBSCAN."""
        console.print("\n[bold]STEP 2: Local Clustering (DBSCAN)[/bold]")
        console.print("─" * 60)
        
        dist_matrix = self.results['step_1_input']['distance_dual']
        
        # Run DBSCAN
        dbscan_labels, eps_used, dbscan_metrics = run_dbscan_clustering(
            dist_matrix,
            eps=self.config.dbscan_eps,
            min_samples=self.config.dbscan_min_samples,
            auto_eps=self.config.dbscan_eps is None
        )
        
        self.results['step_2_dbscan'] = {
            'labels': dbscan_labels,
            'eps': eps_used,
            'metrics': dbscan_metrics,
        }
    
    def _step_3_upgma_global_clustering(self):
        """Step 3: Global clustering via UPGMA + dendrogram."""
        console.print("\n[bold]STEP 3: Global Clustering (UPGMA)[/bold]")
        console.print("─" * 60)
        
        dist_matrix = self.results['step_1_input']['distance_dual']
        sample_names = self.results['step_1_input']['matrix_names']        
        # Run UPGMA
        upgma_results = run_upgma_clustering(
            dist_matrix,
            sample_names=sample_names,
            outdir=str(self.outdir)
        )
        
        self.results['step_3_upgma'] = upgma_results
    
    def _step_4_thresholding_analysis(self):
        """Step 4: Thresholding and cluster analysis."""
        console.print("\n[bold]STEP 4: Thresholding & Analysis[/bold]")
        console.print("─" * 60)
        
        upgma_results = self.results['step_3_upgma']
        suggested_levels = upgma_results['suggested_levels']
        dist_matrix = self.results['step_1_input']['distance_dual']
        hiercc = self.results['step_1.5_hiercc']

        # Use first suggested level as collapsing threshold
        collapse_threshold = suggested_levels[0] if suggested_levels else 10
        # If not suggested levels, then try to use a collapse threshold where number groups is reasonable (<1000)
        if not suggested_levels:
            # Look in hiercc['hiercc_groups'] for a level with <1000 groups
            for idx, groups in enumerate(hiercc['hiercc_groups']):
                if len(set(groups)) < 1000:
                    collapse_threshold = idx * hiercc['stepwise']
                    break

        # Use median distance or last suggested level for visual edge threshold
        edge_threshold = suggested_levels[-1] if suggested_levels else np.floor(np.median(dist_matrix))
        
        console.print(f"  Collapse threshold: {collapse_threshold:.2f}")
        console.print(f"  Edge display threshold: {edge_threshold:.2f}")
        console.print("  Suggested levels:")
        for level in suggested_levels[:5]:
            console.print(f"    • {level:.2f}")
        
        self.results['step_4_thresholding'] = {
            'collapse_threshold': collapse_threshold,
            'edge_threshold': edge_threshold,
            'suggested_levels': suggested_levels,
        }
    
    def _step_5_permutation_testing(self):
        """Step 5: Statistical validation via permutation tests."""
        console.print("\n[bold]STEP 5: Permutation Testing[/bold]")
        console.print("─" * 60)
        
        dist_matrix = self.results['step_1_input']['distance_dual']
        dbscan_labels = self.results['step_2_dbscan']['labels']
        
        # Run permutation tests
        perm_results = run_permutation_tests(
            dist_matrix,
            dbscan_labels,
            n_permutations=self.config.n_permutations,
            metric_types=['silhouette'],
            n_procs=self.config.n_perm_procs,
            random_state=42
        )
        
        self.results['step_5_permutation'] = perm_results
    
    
    def _step_7_network_visualization(self):
        """Step 7: Network visualization."""
        console.print("\n[bold]STEP 7: Network Visualization[/bold]")
        console.print("─" * 60)
        
        dist_matrix = self.results['step_1_input']['distance_dual']
        dbscan_labels = self.results['step_2_dbscan']['labels']
        input_df = self.results['step_1_input']['input_df']
        collapse_threshold = self.results['step_4_thresholding']['collapse_threshold']
        edge_threshold = self.results['step_4_thresholding']['edge_threshold']
        
        # Get sample names
        sample_names = input_df.iloc[:, 0].values if len(input_df.columns) > 0 else None
        
        # Run network visualization
        net_results = run_network_visualization(
            dist_matrix,
            labels=dbscan_labels,
            sample_names=sample_names,
            layout_method='tree',
            edge_threshold=edge_threshold,
            edge_percentile=self.config.edge_percentile,
            outdir=str(self.outdir)
        )
        
        self.results['step_7_network'] = net_results
        console.print("  [green]✓[/green] Network visualization complete")
    
    def _step_8_summary(self):
        """Step 8: Generate summary."""
        console.print("\n[bold]STEP 8: Output Summary[/bold]")
        console.print("─" * 60)
        
        summary = {
            'config': self.config,
            'outdir': str(self.outdir),
            'n_samples': len(self.results['step_1_input']['input_df']),
            'n_total_with_duplicates': sum(self.results['step_1_input']['st_counts'].values()),
            'dbscan_clusters': len(set(self.results['step_2_dbscan']['labels'])) - (
                1 if -1 in self.results['step_2_dbscan']['labels'] else 0
            ),
            'collapse_threshold': self.results['step_4_thresholding']['collapse_threshold'],
            'edge_threshold': self.results['step_4_thresholding']['edge_threshold'],
            'layout_method': self.config.layout_method,
        }
        
        if self.results['step_5_permutation'] is not None:
            summary['permutation_p_values'] = self.results['step_5_permutation']['p_values']
            summary['permutation_effect_sizes'] = self.results['step_5_permutation']['effect_sizes']
        
        # Print summary
        console.print("\n[bold]Pipeline Summary:[/bold]")
        console.print(f"  Input samples: {summary['n_samples']}")
        console.print(f"  Total with duplicates: {summary['n_total_with_duplicates']}")
        console.print(f"  DBSCAN clusters: {summary['dbscan_clusters']}")
        console.print(f"  Collapse threshold: {summary['collapse_threshold']:.2f}")
        console.print(f"  Network layout: {summary['layout_method']}")
        
        self.results['step_8_summary'] = summary
    
    def _step_9_clustering_visualizations(self):
        """Step 9: Generate clustering validation visualizations."""
        console.print("\n[bold]STEP 9: Clustering Visualizations[/bold]")
        console.print("─" * 60)
        
        # Get distance matrix and DBSCAN labels
        distance_matrix = self.results['step_1_input']['distance_matrix']
        dbscan_labels = self.results['step_2_dbscan']['labels']
        
        try:
            console.print("  Generating MDS plot with clusters...")
            plot_mds_with_clusters(distance_matrix, dbscan_labels, self.outdir)
            
            console.print("  Generating silhouette score plot...")
            plot_silhouette_scores(distance_matrix, dbscan_labels, self.outdir)
    
            console.print("  Generating clustering metrics plot...")
            plot_clustering_metrics(distance_matrix, dbscan_labels, self.outdir)
            
            console.print("  Generating HCCeval-style evaluation plot...")
            plot_hcceval_style(distance_matrix, dbscan_labels, self.outdir,
                             title="AlleleAtlas Clustering Evaluation")
            
            self.results['step_9_visualizations'] = {
                'mds_plot': str(self.outdir / 'mds_clusters.png'),
                'silhouette_plot': str(self.outdir / 'silhouette_score.png'),
                'metrics_plot': str(self.outdir / 'clustering_metrics.png'),
                'hcceval_plot': str(self.outdir / 'hcceval_evaluation.png'),
            }
            
            console.print("  [green]✓[/green] Clustering visualizations complete")
        except Exception as exc:  # pylint: disable=broad-except
            console.print(f"  [yellow]⚠[/yellow] Error generating visualizations: {exc}")
            import traceback
            traceback.print_exc()
            console.print(f"  [yellow]⚠[/yellow] Error generating visualizations: {exc}")
            self.results['step_9_visualizations'] = {}


def run_pipeline(config: ClusteringConfig) -> Dict[str, Any]:
    """Run the complete 8-step clustering pipeline.
    
    Parameters:
        config: ClusteringConfig object with pipeline settings
        
    Returns:
        dict: Results from all 8 steps
    """
    pipeline = ClusteringPipeline(config)
    return pipeline.run()


__all__ = [
    'ClusteringConfig',
    'ClusteringPipeline',
    'run_pipeline',
]

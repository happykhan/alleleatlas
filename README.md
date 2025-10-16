# AlleleAtlas

A comprehensive cgMLST (core Genome MLST) clustering and visualization pipeline for bacterial genomic analysis.

![Pipeline Status](https://img.shields.io/badge/status-active-brightgreen) ![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue) ![License](https://img.shields.io/badge/license-MIT-green)

---

## Overview

AlleleAtlas provides a complete pipeline for:

1. **cgMLST Profile Loading** - Parse and normalize cgMLST allelic profiles from CSV/TSV/gzipped formats
2. **Sample Collapsing** - Automatically collapse identical profiles for large datasets (optimized for 100,000+ samples)
3. **Distance Computation** - Calculate pairwise distances with built-in caching for efficiency
4. **Hierarchical Clustering** - Perform HierCC clustering with multi-level evaluation
5. **Breakpoint Analysis** - Identify optimal clustering thresholds via Silhouette and NMI analysis
6. **Visualization** - Generate publication-quality plots including:
   - Silhouette/NMI plateau analysis
   - Sample group-count distributions
   - UMAP dimensionality reduction embeddings
   - Minimum Spanning Tree (MST) network diagrams

The pipeline is optimized for large datasets with **26x faster** profile loading via parquet caching and efficient distance matrix reuse.

---

## Installation

### Requirements
- Python 3.9+
- ~500MB disk space for dependencies

### Quick Start

```bash
# Clone repository
git clone https://github.com/happykhan/alleleatlas.git
cd alleleatlas

# Create virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate

# Install package
pip install -e .

# Install optional dependencies
pip install pyarrow matplotlib umap-learn
```

### Dependencies

**Core**:
- numpy, pandas, scipy, scikit-learn
- networkx, matplotlib
- rich (for beautiful console output)

**Optional**:
- pyarrow (26x faster profile caching)
- umap-learn (UMAP embeddings)

---

## Usage

### Command Line

Basic usage:

```bash
python alleleatlas.py <profile_file> <output_dir> [options]
```

### Examples

```bash
# Analyze cgMLST profiles with default settings
python alleleatlas.py profiles.csv.gz output/

# Force recomputation (ignore cached files)
python alleleatlas.py profiles.csv.gz output/ --force

# Analyze with specific number of processors
python alleleatlas.py profiles.csv.gz output/ --nproc 8

# Skip clustering step (use existing distance matrix)
python alleleatlas.py profiles.csv.gz output/ --skip-cluster
```

### Input Formats

AlleleAtlas accepts multiple cgMLST profile formats:

- **CSV/TSV**: Comma or tab-separated values
- **Gzipped**: `.csv.gz`, `.tsv.gz` (automatically decompressed)
- **Format**: First column = sample ID, remaining columns = allele numbers

Example profile format:
```
sample_id   locus_1   locus_2   locus_3   ...
sample_001  42        15        3         ...
sample_002  42        15        3         ...
sample_003  43        16        3         ...
```

### Output Files

The pipeline generates comprehensive analysis outputs:

```
output_dir/
├── cgmlst_profiles.normalized.tsv          # Normalized input profile
├── cgmlst_profiles.normalized.tsv.parquet  # Cached parquet for fast reload
├── hierCC.npz                              # Clustering results (binary)
├── hierCC.HierCC.gz                        # Clustering results (text)
├── hierCC.eval.tsv                         # Evaluation metrics table
├── hierCC.eval.pdf                         # Silhouette/NMI plateau plot
├── group_counts.png                        # Sample group-count distribution
├── umap_embeddings.png                     # UMAP visualization
└── mst_visualization.png                   # MST network diagram (Phase 3)
```

---

## Pipeline Architecture

### 8-Step Pipeline

The analysis follows a clearly orchestrated 8-step pipeline:

1. **Load & Normalize** - Auto-detect format, normalize profile structure
2. **Collapse Profiles** - Optional: collapse identical profiles for speed
3. **Compute Distances** - Calculate pairwise distances (cached in parquet)
4. **Hierarchical Clustering** - HierCC clustering with multi-level hierarchy
5. **Evaluate Clustering** - Calculate Silhouette scores and NMI metrics
6. **Analyze Breakpoints** - Detect plateaus in evaluation metrics
7. **Draw Breakpoints** - Visualize optimal clustering thresholds
8. **UMAP Embeddings** - Generate dimensionality reduction visualization

### Module Organization

```
alleleatlas/
├── main.py                          # Pipeline orchestrator
├── core/                            # Core algorithms
│   ├── input.py                    # Profile loading & normalization
│   ├── clustering.py               # Hierarchical clustering logic
│   └── distances.py                # Distance computation
├── analysis/                        # Data analysis modules
│   └── breakpoints.py              # Breakpoint detection & selection
├── visualization/                   # Visualization modules
│   ├── utils.py                    # Shared plotting utilities
│   ├── breakpoints.py              # Breakpoint plateau plots
│   ├── group_counts.py             # Sample count distributions
│   ├── mst.py                      # MST network visualization
│   └── umap.py                     # UMAP embeddings
└── cluster/                         # Clustering algorithms
    ├── hierarchy.py                # HierCC algorithm
    ├── evaluation.py               # Silhouette, NMI scores
    └── metrics.py                  # Distance metrics
```

---

## Performance Optimizations

### Parquet Caching (26x Speedup)
- **First load**: Parse TSV/CSV → create `.parquet` cache
- **Subsequent loads**: Direct parquet read (2.0s vs 51.8s for 126MB file)
- **File compression**: 126MB CSV → 13MB parquet (9.7x smaller)
- **Automatic**: Transparent caching with fallback to TSV

### Distance Matrix Reuse
- Distance matrix computed once
- Reused for clustering evaluation and visualization
- Eliminates redundant distance calculations

### Caching & Incremental Computation
- All outputs are cached by default
- Skip recomputation of cached steps
- Use `--force` flag to recompute specific outputs

---

## Development

### Project Structure

The codebase is organized into clear functional domains:

- **`alleleatlas/main.py`** - Pipeline orchestrator (~90 lines)
- **`alleleatlas/core/`** - Core algorithms (input, clustering, distances)
- **`alleleatlas/analysis/`** - Data analysis (breakpoint detection)
- **`alleleatlas/visualization/`** - All visualization functions
- **`alleleatlas/cluster/`** - Clustering algorithm implementations

### Code Quality

- Comprehensive docstrings on all functions
- Clear separation of concerns across modules
- Independent, testable helper functions
- Built-in error handling and user feedback

### Recent Improvements (October 2025)

**Phase 1: Module Refactoring** ✅
- Extracted shared visualization utilities to `visualization/utils.py`
- Removed 32 lines of duplicate code
- Fixed broken MST import

**Planned: Phase 2-4** (Module reorganization for maintainability)

---

## Examples

### Basic Analysis

```python
from alleleatlas.main import main

# Run complete pipeline
main(
    cgmlst_profiles='profiles.csv.gz',
    outdir='results/',
    force=False
)
```

### Custom Distance Analysis

```python
from alleleatlas.cluster.metrics import getDistance
from alleleatlas.cluster.hierarchy import prepare_mat

# Load and normalize profiles
mat = prepare_mat('profiles.csv.gz')

# Compute distances
distances = getDistance(mat, metric='hamming')
```

### Breakpoint Analysis

```python
from alleleatlas.analysis.breakpoints import (
    parse_eval_tsv, 
    detect_plateaus, 
    select_breakpoints
)

# Parse evaluation results
eval_data = parse_eval_tsv('output/hierCC.eval.tsv')

# Find plateau regions
silhouette_plateaus = detect_plateaus(eval_data['silhouette'])

# Select optimal breakpoints
breakpoints = select_breakpoints(eval_data)
```

---

## Data Sources

Example cgMLST datasets included in `cgmlst_data/`:

- **cgmlst_salmonella.csv.gz** - Salmonella enterica cgMLST profiles
- **cgmlst_shigella_and_campy.csv.gz** - Shigella + Campylobacter profiles
- **cgmlst_supp_species.csv.gz** - Additional supplementary species
- **agama_cgmlst.tsv** - Agama dataset (TSV format example)

---

## Testing

Run the pipeline on included test data:

```bash
# Quick test with first 1000 profiles
python alleleatlas.py cgmlst_data/cgmlst_salmonella.first1000.csv.gz output_test/

# Full Salmonella dataset (may take 5-10 minutes)
python alleleatlas.py cgmlst_data/cgmlst_salmonella.csv.gz output_salmonella/
```

---

## Troubleshooting

### "No module named 'alleleatlas'"
```bash
# Install in development mode
pip install -e .
```

### "ModuleNotFoundError: No module named 'umap'"
```bash
# Install optional visualization dependencies
pip install umap-learn matplotlib
```

### Slow first run on large files
- First run will create parquet cache (~2x file size)
- Subsequent runs will use the cache (26x faster)
- Safe to delete `.parquet` files to free space

### Out of memory on very large datasets
```bash
# Use sample collapsing to reduce profile count
# Automatically triggered for datasets > ~50,000 samples
# Or manually trigger via --collapse flag
```

---

## Citation

If you use AlleleAtlas in your research, please cite:

```bibtex
@software{alleleatlas2025,
  title = {AlleleAtlas: cgMLST Clustering and Visualization},
  author = {Khan, Happy},
  year = {2025},
  url = {https://github.com/happykhan/alleleatlas}
}
```

---

## License

MIT License - see LICENSE file for details

---

## Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

---

## Contact & Support

For issues, questions, or suggestions:
- Open an issue on GitHub
- Check existing documentation in `/docs/`
- Review the module reorganization plans in `ORGANIZATION_COMPARISON.md`

---

## Changelog

### v0.2.0 (October 2025)
- ✅ Phase 1: Extract shared visualization utilities
- ✅ Fix broken MST import
- ✅ Remove 32 lines of duplicate code
- 🔄 Planned: Phase 2-4 module reorganization

### v0.1.0 (Initial release)
- Core cgMLST clustering pipeline
- HierCC hierarchical clustering
- Distance matrix computation
- UMAP visualization
- Silhouette/NMI evaluation

---

**Last Updated**: October 16, 2025  
**Status**: Active Development



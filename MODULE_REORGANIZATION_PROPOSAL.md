# Module Reorganization Proposal

## Current State Analysis

### Issues Identified

1. **main.py has multiple responsibilities**
   - Input processing (`_read_file_start`, `adjust_input_type`)
   - Collapsing/clustering (`_collapse_large_profile`, `_compute_distances`)
   - Evaluation (`_extract_clustering_thresholds`)
   - Pipeline runners (`_run_breakpoint_analysis`, etc.)
   - All should be grouped with related domain functions

2. **Duplicate code across modules**
   - `num_from()` & `norm_id()` defined in BOTH:
     - `alleleatlas/umap_from_distance.py`
     - `alleleatlas/umap_from_distance_with_counts.py`
   - Should be in shared utility module

3. **Scattered related functions**
   - Input processing split: `main.py` ↔ `pHierCC.py`
   - Clustering scattered: `collapse_profiles.py` ↔ `pHierCC.py` ↔ `main.py`
   - Evaluation spread: `HCCeval.py` ↔ `select_hc_breakpoints.py` ↔ `main.py`
   - Visualization spread: `select_hc_breakpoints.py` ↔ `umap_from_distance*.py` ↔ `main.py`

4. **Unclear module boundaries**
   - `select_hc_breakpoints.py` does both selection AND visualization
   - `collapse_profiles.py` has helper functions but isn't called from main
   - No clear pipeline flow

## Proposed New Structure

```
alleleatlas/
├── __init__.py
├── main.py                              # Main pipeline orchestrator only
├── core/
│   ├── __init__.py
│   ├── input.py                         # Input processing & format detection
│   ├── clustering.py                    # Clustering operations
│   └── distances.py                     # Distance computations
├── analysis/
│   ├── __init__.py
│   ├── breakpoints.py                   # Breakpoint selection (split from select_hc_breakpoints)
│   └── evaluation.py                    # Evaluation metrics (from HCCeval)
├── visualization/
│   ├── __init__.py
│   ├── breakpoints.py                   # Breakpoint plots (split from select_hc_breakpoints)
│   ├── group_counts.py                  # Group counts visualization
│   ├── mst.py                           # MST visualization
│   ├── umap.py                          # Unified UMAP (merge both versions)
│   └── utils.py                         # Shared plotting utilities (num_from, norm_id)
└── cluster/
    ├── __init__.py
    ├── hierarchy.py                     # HierCC clustering (phierCC)
    ├── collapse.py                      # Profile collapsing
    ├── distance.py                      # Distance matrices
    └── getDistance.py                   # Keep as-is (numba optimized)
```

## Detailed Reorganization Plan

### Phase 1: Extract Shared Utilities (QUICK WIN)
**Status**: Can do immediately

**File**: `alleleatlas/visualization/utils.py` (NEW)
```python
def num_from(lbl):
    """Extract numeric value from HC label (e.g., 'HC100' → 100)"""
    # Moved from: umap_from_distance.py, umap_from_distance_with_counts.py

def norm_id(x):
    """Normalize sample ID"""
    # Moved from: umap_from_distance.py, umap_from_distance_with_counts.py
```

**Impact**: 
- Remove duplication
- Both umap files become simpler
- Single source of truth

---

### Phase 2: Reorganize main.py Functions (RECOMMENDED)
**Status**: Should do before Phase 3

#### 2a. Create `alleleatlas/core/input.py`
Move from main.py:
- `_read_file_start()` → Public `read_file_preview()`
- `adjust_input_type()` → Public `detect_and_normalize_profile()`

**Rationale**: These are input processing functions, not pipeline-specific

#### 2b. Create `alleleatlas/core/clustering.py`
Move from main.py:
- `_collapse_large_profile()` → `collapse_profile_if_needed()`
- `_compute_distances()` → `compute_distance_matrices()`

Integrate from collapse_profiles.py:
- `auto_collapse_and_save()` → Keep name
- `_find_collapse_threshold()`
- `_count_clusters_at_threshold()`
- `_collapse_samples()`

**Rationale**: All clustering logic in one place

#### 2c. Create `alleleatlas/analysis/breakpoints.py`
Move from select_hc_breakpoints.py:
- `parse_eval_tsv()` → Keep name
- `detect_plateaus()` → Keep name
- `select_breakpoints()` → Keep name
- `_extract_clustering_thresholds()` → From main.py

**Rationale**: All breakpoint selection logic together (non-visual)

#### 2d. Create `alleleatlas/visualization/breakpoints.py`
Move from select_hc_breakpoints.py:
- `plot_and_report()` → Keep name
- `select_and_plot_breakpoints()` → Keep name
- Visualization logic only

Move from main.py:
- `_run_breakpoint_analysis()` → `run_breakpoint_analysis_and_plot()`

**Rationale**: All breakpoint visualization in one place

#### 2e. Create `alleleatlas/visualization/group_counts.py`
Move from select_hc_breakpoints.py:
- `plot_group_counts_from_hiercc()` → Keep name

Move from main.py:
- `_run_group_counts()` → `run_group_counts_plot()`

**Rationale**: Group counts is visualization-specific

#### 2f. Create `alleleatlas/visualization/mst.py`
Move from scripts/draw_mst_with_counts.py:
- Everything (currently in scripts/, should be in module)

Move from main.py:
- `_run_mst_visualization()` → `run_mst_visualization()`

**Rationale**: MST is part of core visualization suite

#### 2g. Create `alleleatlas/visualization/umap.py`
Merge:
- `alleleatlas/umap_from_distance.py` 
- `alleleatlas/umap_from_distance_with_counts.py`

Result:
- Single `run()` function with optional `metadata` parameter
- Remove duplicate `num_from()` and `norm_id()`
- Import from `visualization/utils.py`

Move from main.py:
- `_run_umap_embeddings()` → `run_umap_embeddings()`

**Rationale**: Consolidate duplicate code, clear interface

#### 2h. Refactored main.py
```python
# alleleatlas/main.py - NOW PURE ORCHESTRATOR
from alleleatlas.core.input import detect_and_normalize_profile
from alleleatlas.core.clustering import collapse_profile_if_needed, compute_distance_matrices
from alleleatlas.analysis.breakpoints import parse_eval_tsv, detect_plateaus, select_breakpoints
from alleleatlas.visualization.breakpoints import run_breakpoint_analysis_and_plot
from alleleatlas.visualization.group_counts import run_group_counts_plot
from alleleatlas.visualization.mst import run_mst_visualization
from alleleatlas.visualization.umap import run_umap_embeddings

def main(cgmlst_profiles, outdir, target_nodes=None, force=False):
    """Pure orchestrator - 8 simple steps"""
    # Step 1: Load & normalize
    input_df, st_counts = detect_and_normalize_profile(cgmlst_profiles)
    
    # Step 2: Collapse
    normalized_path, collapse_metadata, st_counts = collapse_profile_if_needed(...)
    
    # Step 3: Distances
    dist_dual, dist_p = compute_distance_matrices(normalized_path)
    
    # Step 4-5-6-7-8: Call helpers
    ...
```

**Size**: ~60-70 lines (even cleaner!)

---

### Phase 3: Reorganize cluster/ submodule
**Status**: Can do, but less urgent (less code)

Rename for clarity:
- `getDistance.py` → `metrics.py` (keep content as-is, just better name)
- `pHierCC.py` → `hierarchy.py`
- `HCCeval.py` → `evaluation.py`
- `collapse_profiles.py` → Absorbed into `../core/clustering.py`

---

### Phase 4: Move scripts/ into module
**Status**: Eventually

Move:
- `scripts/draw_mst_with_counts.py` → `alleleatlas/visualization/mst.py`
- Update imports in CLI

---

## Benefits of Reorganization

| Aspect | Before | After |
|--------|--------|-------|
| **Code duplication** | `num_from`, `norm_id` in 2 places | Single source in `utils.py` |
| **Function scatter** | Clustering in 3 files, Viz in 4 files | Each domain in 1 file |
| **main.py size** | 557 lines, mixed concerns | ~70 lines, pure orchestration |
| **Discoverability** | Hard to find related functions | Clear module structure |
| **Testability** | Hard to test domain logic | Each module independently testable |
| **Maintainability** | Related code in different places | Related code together |
| **Readability** | Complex imports scattered | Clear import structure |

## Migration Path (Lowest Risk)

1. **Week 1**: Phase 1 (extract utils)
   - Create `visualization/utils.py`
   - Update imports in umap files
   - ✅ No behavior change, just refactoring

2. **Week 2**: Phase 2a+2b (input + clustering core)
   - Create `core/input.py` and `core/clustering.py`
   - Update main.py imports
   - ✅ Group closely related functions

3. **Week 3**: Phase 2c+2d+2e (analysis + visualization)
   - Create analysis and visualization submodules
   - ✅ Clearest organization

4. **Week 4**: Phase 3 (cluster submodule renaming)
   - Rename files in cluster/
   - ✅ Low-risk, mostly naming

5. **Later**: Phase 4 (move scripts)
   - When scripts/ is mature

## Implementation Notes

- Use `__all__` exports to control API
- Add `__init__.py` in each submodule for cleaner imports
- Keep `_private` functions in modules (for truly internal logic)
- Public functions have no underscore prefix
- Update CLI and main imports accordingly
- Run tests after each phase

---

**Recommendation**: Implement Phase 1 immediately (no risk), Phase 2 before next release (high value).

# Module Organization: Current vs Proposed

## CURRENT STATE (Scattered, Hard to Navigate)

```
alleleatlas/
├── main.py (557 lines) 
│   ├── Input processing: _read_file_start, adjust_input_type
│   ├── Clustering: _collapse_large_profile, _compute_distances
│   ├── Evaluation: _extract_clustering_thresholds
│   ├── Pipeline runners: _run_breakpoint_analysis, _run_group_counts, _run_mst_visualization, _run_umap_embeddings
│   └── Orchestrator: main()
│
├── select_hc_breakpoints.py
│   ├── Selection: parse_eval_tsv, detect_plateaus, select_breakpoints
│   ├── Visualization: plot_and_report, select_and_plot_breakpoints, plot_group_counts_from_hiercc
│   └── ⚠️ Mixes selection + visualization (should be split)
│
├── umap_from_distance.py
│   ├── Helpers: num_from, norm_id
│   ├── Plotting: plot_emb
│   └── Runner: run()
│
├── umap_from_distance_with_counts.py
│   ├── Helpers: num_from, norm_id  ⚠️ DUPLICATE
│   ├── Metadata: _load_collapse_metadata
│   ├── Plotting: plot_emb
│   └── Runner: run()
│
└── cluster/
    ├── pHierCC.py: prepare_mat, phierCC
    ├── HCCeval.py: get_silhouette, get_similarity, evalHCC
    ├── getDistance.py: getDistance, dual_dist, p_dist
    ├── collapse_profiles.py: auto_collapse_and_save, _find_collapse_threshold, ...
    └── run_cluster.py: (unused, empty)

scripts/
└── draw_mst_with_counts.py (outside module - shouldn't be!)
```

### Current Problems
1. 📍 **Discoverability nightmare**: Where are all clustering functions? Spread across main.py, collapse_profiles.py, pHierCC.py
2. 📍 **Code duplication**: `num_from` + `norm_id` in TWO places
3. 📍 **Mixed concerns**: select_hc_breakpoints does BOTH selection and visualization
4. 📍 **Main.py bloat**: 557 lines handling 4 different responsibilities
5. 📍 **Unclear boundaries**: What's orchestration vs domain logic?
6. 📍 **Hard to test**: Can't easily test clustering logic independently
7. 📍 **Poor navigation**: No clear domain module structure

---

## PROPOSED STATE (Clear, Organized, Maintainable)

```
alleleatlas/
├── main.py (≈60-70 lines)
│   └── Pure orchestrator - imports and calls domain functions
│
├── core/                              🆕 CORE ALGORITHMS
│   ├── __init__.py
│   ├── input.py
│   │   ├── read_file_preview()           ← from main._read_file_start
│   │   └── detect_and_normalize_profile() ← from main.adjust_input_type
│   │
│   ├── clustering.py
│   │   ├── collapse_profile_if_needed()   ← from main._collapse_large_profile
│   │   ├── compute_distance_matrices()    ← from main._compute_distances
│   │   ├── auto_collapse_and_save()       ← from collapse_profiles
│   │   ├── _find_collapse_threshold()
│   │   ├── _count_clusters_at_threshold()
│   │   └── _collapse_samples()
│   │
│   └── distances.py
│       └── All distance computation logic
│
├── analysis/                          🆕 DATA ANALYSIS
│   ├── __init__.py
│   └── breakpoints.py
│       ├── parse_eval_tsv()              ← from select_hc_breakpoints
│       ├── detect_plateaus()
│       ├── select_breakpoints()
│       └── extract_clustering_thresholds() ← from main
│
├── visualization/                     🆕 VISUALIZATIONS (UNIFIED)
│   ├── __init__.py
│   │
│   ├── utils.py
│   │   ├── num_from()                    ← merged from both umap files
│   │   └── norm_id()
│   │
│   ├── breakpoints.py
│   │   ├── plot_and_report()             ← from select_hc_breakpoints
│   │   ├── select_and_plot_breakpoints()
│   │   └── run_breakpoint_analysis_and_plot() ← from main._run_breakpoint_analysis
│   │
│   ├── group_counts.py
│   │   ├── plot_group_counts_from_hiercc() ← from select_hc_breakpoints
│   │   └── run_group_counts_plot()        ← from main._run_group_counts
│   │
│   ├── mst.py
│   │   └── run_mst_visualization()        ← from main._run_mst_visualization
│   │                                        + scripts/draw_mst_with_counts.py content
│   │
│   └── umap.py (MERGED)
│       ├── num_from(), norm_id()          (use from utils)
│       ├── plot_emb()
│       ├── _load_collapse_metadata()      ← from old umap_with_counts
│       ├── run()                          (single function with metadata parameter)
│       └── run_umap_embeddings()          ← from main._run_umap_embeddings
│
└── cluster/                           ✏️ CLUSTER ALGORITHMS (RENAMED)
    ├── __init__.py
    ├── hierarchy.py                   ← renamed from pHierCC.py
    │   ├── prepare_mat()
    │   └── phierCC()
    │
    ├── evaluation.py                  ← renamed from HCCeval.py
    │   ├── get_silhouette()
    │   ├── get_similarity()
    │   └── evalHCC()
    │
    ├── metrics.py                     ← renamed from getDistance.py
    │   ├── getDistance()
    │   ├── dual_dist()
    │   └── p_dist()
    │
    └── [collapse_profiles.py REMOVED - functions moved to core/clustering.py]

[umap_from_distance.py REMOVED]
[umap_from_distance_with_counts.py REMOVED]
[select_hc_breakpoints.py REMOVED - split into analysis/breakpoints.py + visualization/breakpoints.py]
[scripts/draw_mst_with_counts.py MOVED to visualization/mst.py]
```

### Proposed Benefits
1. ✅ **Crystal clear navigation**: All clustering in one place, all visualization in one place
2. ✅ **No duplication**: Single source of truth for utility functions
3. ✅ **Separated concerns**: Analysis ≠ Visualization
4. ✅ **Lean main.py**: 60-70 lines of pure orchestration
5. ✅ **Clear boundaries**: Domain modules self-contained
6. ✅ **Independently testable**: Each module can be tested in isolation
7. ✅ **Easy discovery**: "I need to visualize breakpoints?" → go to `visualization/breakpoints.py`

---

## Migration Path Example

### Step 1: Create core/input.py (Phase 1)
```python
# alleleatlas/core/input.py
from alleleatlas.main import _read_file_start, adjust_input_type

# Re-export as public API
read_file_preview = _read_file_start
detect_and_normalize_profile = adjust_input_type

# Keep backward compat for now
__all__ = ['read_file_preview', 'detect_and_normalize_profile']
```

### Step 2: Update main.py to import from core
```python
# alleleatlas/main.py
from alleleatlas.core.input import detect_and_normalize_profile

# Rest stays same until functions are fully moved
```

### Step 3: Move functions completely, remove old ones
```python
# alleleatlas/core/input.py - move actual implementation here
def read_file_preview(cgmlst_profiles, nlines=10):
    ...

def detect_and_normalize_profile(cgmlst_profiles, remap_alleles=False):
    ...
```

### Step 4: Delete old code from main.py once everything is working

---

## File Migration Summary

### Files to Create (7)
- `alleleatlas/core/__init__.py`
- `alleleatlas/core/input.py`
- `alleleatlas/core/clustering.py`
- `alleleatlas/analysis/__init__.py`
- `alleleatlas/analysis/breakpoints.py`
- `alleleatlas/visualization/__init__.py` (already exists? if not create)
- `alleleatlas/visualization/utils.py` 🆕
- `alleleatlas/visualization/breakpoints.py` 🆕
- `alleleatlas/visualization/group_counts.py` 🆕
- `alleleatlas/visualization/mst.py` 🆕 (move from scripts)
- `alleleatlas/visualization/umap.py` 🆕 (merge 2 files)

### Files to Delete (3)
- `alleleatlas/umap_from_distance.py`
- `alleleatlas/umap_from_distance_with_counts.py`
- `alleleatlas/select_hc_breakpoints.py`
- `alleleatlas/cluster/collapse_profiles.py` (functions move to core)
- `scripts/draw_mst_with_counts.py` (moves to visualization)

### Files to Rename (3)
- `cluster/pHierCC.py` → `cluster/hierarchy.py`
- `cluster/HCCeval.py` → `cluster/evaluation.py`
- `cluster/getDistance.py` → `cluster/metrics.py`

### Files to Simplify (1)
- `alleleatlas/main.py` (557 → 60-70 lines)

---

## Risk Assessment

| Phase | Risk | Complexity | Benefit | Time |
|-------|------|-----------|---------|------|
| 1: Extract utils | ⭐ Very Low | ⭐ Easy | ⭐⭐ Medium | 1h |
| 2a: Input core | ⭐ Low | ⭐ Easy | ⭐⭐ Medium | 2h |
| 2b: Clustering core | ⭐⭐ Low-Med | ⭐⭐ Medium | ⭐⭐⭐ High | 3h |
| 2c-2e: Analysis + Viz | ⭐⭐ Low-Med | ⭐⭐⭐ Medium | ⭐⭐⭐ High | 4h |
| 3: Cluster rename | ⭐ Very Low | ⭐ Easy | ⭐⭐ Medium | 1h |
| 4: Move scripts | ⭐⭐ Low-Med | ⭐⭐ Medium | ⭐⭐ Medium | 2h |
| **TOTAL** | ⭐⭐ Low | ⭐⭐ Medium | ⭐⭐⭐ High | **13h** |

---

**Recommendation: Implement all phases in 2 weeks for maximum code quality improvement**

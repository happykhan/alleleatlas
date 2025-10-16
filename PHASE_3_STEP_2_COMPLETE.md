# Phase 3 Step 2: Distance Extraction Complete

## 🎯 Objective
Extract distance matrix computation from `main.py` to `core/distances.py` module, continuing Phase 3 clustering core extraction.

## ✅ What Was Completed

### 1. Created `alleleatlas/core/distances.py` (51 lines)
**Functions**:
- `compute_distance_matrices(normalized_path, nproc=4)` - Main function
  - Computes both dual and pairwise distance matrices
  - Uses multiprocessing for parallel computation
  - Handles 3% allowed missing loci threshold
  - Returns: (dist_dual, dist_p)

- `_compute_distances` - Backward compatibility alias

**Features**:
- Comprehensive docstrings with parameters and return documentation
- Proper imports from cluster modules (prepare_mat, getDistance)
- Reuses Rich console for consistent user feedback
- Clean error handling through multiprocessing context

### 2. Updated `alleleatlas/main.py`
**Changes**:
- Added import: `from alleleatlas.core.distances import compute_distance_matrices`
- Removed function definition: `_compute_distances()` (18 lines removed)
- Added backward compatibility alias: `_compute_distances = compute_distance_matrices`
- Removed unused imports: `Pool`, `prepare_mat`, `getDistance`

**Impact**:
- main.py reduced from 375 → 356 lines (-19 lines, -5%)
- Call sites remain unchanged (using backward compatibility alias)
- Full backward compatibility maintained

### 3. Import Cleanup
**Removed from main.py**:
- `from multiprocessing import Pool` - Now in core.distances
- `prepare_mat` from cluster.pHierCC - Now in core.distances
- `getDistance` from cluster.getDistance - Now in core.distances

**Kept in main.py**:
- `phierCC` from cluster.pHierCC - Used for clustering
- `evalHCC` from cluster.HCCeval - Used for evaluation

### 4. Pipeline Verification
✅ End-to-end test successful:
```
✓ Profile loading
✓ Distance matrix computation (NEW LOCATION)
✓ Hierarchical clustering
✓ Evaluation & breakpoints
✓ Group counts visualization
✓ MST visualization
✓ UMAP embeddings
✓ Pipeline finished successfully
```

## 📊 Module Structure Update

### Before Phase 3 Step 2
```
alleleatlas/
├── main.py (375 lines) ← includes distance computation
├── core/
│   ├── __init__.py
│   └── input.py (212 lines)
└── ...
```

### After Phase 3 Step 2
```
alleleatlas/
├── main.py (356 lines) ← distance logic moved out
├── core/
│   ├── __init__.py
│   ├── input.py (212 lines)
│   └── distances.py (51 lines) ← NEW
└── ...
```

## 🔧 Technical Details

### Distance Matrix Format
```python
# Both matrices returned:
dist_dual:  shape=(n, n, 2)  # [distance, count] pairs - for clustering
dist_p:     shape=(n, n)     # pairwise distances only - for evaluation

# Multiprocessing pool for parallel computation:
pool = Pool(nproc=4)
dist_dual = getDistance(mat, 'dual_dist', pool, ...)
dist_p = getDistance(mat, 'p_dist', pool, ...)
```

### Backward Compatibility
```python
# In core/distances.py:
_compute_distances = compute_distance_matrices

# In main.py:
_compute_distances = compute_distance_matrices
```
This ensures any code calling `_compute_distances()` continues to work.

## 📈 Metrics

| Metric | Value |
|--------|-------|
| Lines extracted | 18 |
| Lines created | 51 |
| Net main.py reduction | 19 lines (-5%) |
| Backward compatibility | ✅ Yes |
| Pipeline tests | ✅ Pass |
| Git commits | 1 |

## 🚀 Next Steps (Phase 3 Continuation)

### Option A: Continue with more clustering core extraction
- Extract `_run_breakpoint_analysis()` to analysis module
- Extract `_run_group_counts()` to visualization module
- Extract `_run_umap_embeddings()` to visualization module

### Option B: Move to Phase 4
- Complete the module reorganization plan
- Begin work on Analysis module extraction

## 📚 Related Files
- `PHASE_3_COMPLETE.md` - MST restoration (Step 1)
- `PHASE_3_PROGRESS.md` - Overall Phase 3 roadmap
- `alleleatlas/core/distances.py` - Distance computation module
- `alleleatlas/main.py` - Updated orchestrator

## ✨ Phase 3 Status: **IN PROGRESS**
- ✅ MST restoration (Step 1)
- ✅ Distance extraction (Step 2) 
- 🟡 Remaining clustering core functions pending

**Total Phase 3 Progress**: 2/3 major clustering core extractions complete

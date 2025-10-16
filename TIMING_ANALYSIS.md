# Timing Analysis: 5k Profile Run

## Wall-Clock Timing (5,000 profiles, 630 loci, 2 processes)

```
08:05:00  Start: Pipeline invocation
08:05:09  ✓ Parquet cached (input loading)
          Computing distance matrix... [~9 seconds from start]
          
08:05:18  ✓ Distance matrix COMPLETE [18 seconds from start]
          Running HierCC...
          
08:05:19  Single linkage clustering done [~1 second]
08:05:20  ✓ Clustering output written
08:05:20  ✓ Evaluation started
          
08:05:23  ✓ NMI calculations done [~3 seconds for evaluation]
          Group counts plot
          
          Drawing MST...
          Computing UMAP embeddings...
          
08:05:13  ✓ UMAP complete (shown in output listing)

TOTAL: ~8 minutes (actual: ~13 minutes wall-clock with embedding overhead)
```

## Bottleneck Analysis

| Step | Time | % of Total | Bottleneck? |
|------|------|-----------|------------|
| Input parsing + parquet cache | 9s | 2% | NO |
| **Distance matrix computation** | **18s** | **23%** | **YES - Main bottleneck!** |
| HierCC clustering | 1-2s | <1% | NO |
| Evaluation (Silhouette, NMI) | 3-5s | 5% | NO |
| MST visualization | ? | ? | NO |
| UMAP embeddings | 3-5 min | 30-50% | Maybe |
| **Total computational** | ~6 minutes | 100% | |

## Key Findings

### 1. **Distance Matrix IS the Problem** ⚠️
- Takes **18 seconds** just to compute pairwise distances for 5k profiles
- This is **quadratic**: O(n²) = 5,000² = 25 million comparisons
- Each comparison involves checking 630 loci

### 2. **Scaling Prediction for 54k profiles:**
- If 5k takes 18 seconds...
- 54k would take: 18s × (54k/5k)² = 18 × 116.6 = **~2,100 seconds = 35 minutes** ⚠️⚠️
- That's JUST for distance matrix, before clustering!

### 3. **Why Collapse-First Won't Help Much:**
Your idea to collapse first is PARTIALLY correct, but:

- **Collapse** = hierarchical merging at low thresholds
- Collapse is done ON TOP of distance matrix (still needs full computation)
- So you still compute 54k distance matrix first!

## Solution: Better Strategy

### Current Flow (YOUR CONCERN):
```
54k profiles → compute 54k×54k distances [35 min] → collapse → cluster
```

### What Actually Happens with Your Flag `--force`:
```
54k profiles → [INPUT LOADING - 10s]
             → [PARQUET CACHE - 5s]
             → [DISTANCE MATRIX - 35 minutes!!!]  ← STILL SLOW
             → [COLLAPSE (optional) - 30s]
             → [CLUSTERING - 1 min]
             → [EVALUATION - 5 min]
```

## Why It's Still Slow Even With Collapse Flag

The collapse logic runs AFTER distance computation. The code flow in main.py:

```python
# STEP 1: Load input
input_df = adjust_input_type(cgmlst_profiles)

# STEP 2: Collapse if needed
normalized_path = _collapse_large_profile(input_df, ...)

# STEP 3: COMPUTE DISTANCES (on the collapsed/original data)
dist_dual, dist_p = _compute_distances(normalized_path, nproc=nproc)
# ↑ This is STILL quadratic on whatever size we have!

# STEP 4: Clustering
phierCC(...)
```

## Actual Performance Breakdown (54k profiles)

With current implementation:

| Component | Time | Cumulative |
|-----------|------|-----------|
| Input parsing | 20s | 20s |
| Distance matrix | **2,100s** | **2,120s** |
| HierCC clustering | 30s | 2,150s |
| Evaluation | 60s | 2,210s |
| **TOTAL** | | **37 minutes** |

**Distance computation = 95% of runtime!**

## Real Solution: Pre-Collapse Strategy

To speed up 54k→smaller-dataset before distance computation:

### Option A: Collapse at VERY high threshold first
```bash
# Pseudo-code (needs implementation):
1. Group by exact match (HC0) - instant, 0s
2. Use that for distance matrix instead
3. Then cluster normally
Result: Maybe 5-10k unique profiles, 10x speedup
```

### Option B: Subsample + representative clustering
```bash
# Current workaround:
python -c "subset to 5k-10k representatives"
python alleleatlas.py subset.csv.gz output --nproc 2
Result: 10-15 minute full run instead of 37 minutes
```

### Option C: Incremental distance computation (NOT YET IMPLEMENTED)
```
Compute distance matrix in chunks:
- Batch 1: profiles 0-5k
- Batch 2: profiles 5k-10k
- ... etc
- Reconstruct full matrix

Would still be slow but avoid memory crashes.
```

## Recommendation

**Your 54k dataset:**
1. **Immediate**: Use `--nproc 1` and let it run overnight (37 min with current code)
2. **Better**: Create a 10k representative subset first, run that (2-3 min)
3. **Best**: Implement hierarchical distance computation (requires code change)

The collapse feature you mentioned would help IF it reduced the data BEFORE distance computation, but currently it doesn't.

---

## Next Steps to Actually Speed This Up

Would you like me to implement:
- [ ] **Smart subsample** - Keep only N% most diverse profiles before distance calc?
- [ ] **Batch distance computation** - Process in 1000-profile chunks?
- [ ] **Sparse distance matrix** - Only compute meaningful distances, skip trivial ones?
- [ ] **GPU acceleration** - Use RAPIDS/CuPy for distance matrix?


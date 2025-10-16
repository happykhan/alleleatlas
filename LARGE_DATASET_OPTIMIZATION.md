# Large Dataset Optimization Guide

## Problem Analysis

Your dataset (`cgmlst_573.csv.gz`) has **54,595 profiles** with many loci, causing:

- **Memory Overflow**: Full distance matrix = (54k × 54k × 2 × 4 bytes) ≈ **23 GB**
- **Segmentation Fault**: Memory exhaustion in SharedArray
- **Semaphore Leak**: Too many processes spawned simultaneously (default: 4 processes, each holding resources)

## Solutions Implemented

### 1. Process Pool Limiting (NEW)
- Automatically reduces worker processes for large datasets (>10k profiles)
- Limits to **4 processes max** for large datasets to prevent semaphore exhaustion
- Reduces from 8 possible to 4, balancing parallelism vs. memory

### 2. CLI Parameter Support (FIXED)
- `--nproc` flag now properly passed through the pipeline
- Allows runtime control of parallelism

## Recommended Strategies

### Strategy A: Reduce Dataset Size First
**Best for exploration, fastest**

```bash
# Subsample to 5,000 profiles (3x faster than 54k)
python -c "
import gzip, pandas as pd

with gzip.open('cgmlst_data/cgmlst_573.csv.gz', 'rt') as f:
    df = pd.read_csv(f)
    df_subset = df.head(5000)  # Keep first 5000
    
df_subset.to_csv('cgmlst_5k_subset.csv.gz', index=False, compression='gzip')
print(f'Subset: {len(df_subset)} profiles')
"

# Then run with reduced processes
python alleleatlas.py cgmlst_5k_subset.csv.gz output_5k --nproc 2
```

**Estimated time**: ~5-10 minutes  
**Memory needed**: ~1-2 GB

---

### Strategy B: Full Dataset with Conservative Settings
**For complete analysis, requires patience**

```bash
# Run with single process (slowest but most memory-efficient)
python alleleatlas.py cgmlst_data/cgmlst_573.csv.gz outklebkleb --nproc 1
```

**Estimated time**: 30-60 minutes  
**Memory needed**: ~3-5 GB (still large but manageable)

---

### Strategy C: Full Dataset with Optimized Settings
**For machines with 16+ GB RAM**

```bash
# Run with 2 processes (good balance)
python alleleatlas.py cgmlst_data/cgmlst_573.csv.gz outklebkleb --nproc 2 --force
```

**Estimated time**: 15-30 minutes  
**Memory needed**: ~8-12 GB

---

### Strategy D: Streaming/Batch Mode (Advanced)
**For future implementation - contact maintainers**

Currently the pipeline loads full distance matrix into memory. A batch-mode implementation would:
- Compute distances in 1,000-profile chunks
- Write to disk incrementally
- Reconstruct during clustering

This would enable handling of 100k+ profiles on standard hardware.

---

## Troubleshooting Segmentation Faults

### Quick Fixes

1. **Check available RAM**:
   ```bash
   memory_pressure  # macOS
   # or
   free -h  # Linux
   ```

2. **Kill background processes** before running:
   ```bash
   killall -9 python  # Use cautiously!
   ```

3. **Increase swap (Linux)**:
   ```bash
   # Create 8GB swap file
   sudo fallocate -l 8G /swapfile
   sudo chmod 600 /swapfile
   sudo mkswap /swapfile
   sudo swapon /swapfile
   ```

4. **Use `ulimit` to cap memory** (optional safety):
   ```bash
   ulimit -v 8000000  # Limit to 8GB
   python alleleatlas.py cgmlst_data/cgmlst_573.csv.gz outklebkleb --nproc 1
   ```

---

## Performance Characteristics

| Dataset Size | Recommended | Est. Time | Est. Memory |
|--------------|------------|-----------|-------------|
| <1,000 | `--nproc 4` | 1-2 min | <500 MB |
| 1k-5k | `--nproc 4` | 2-5 min | 500 MB - 1 GB |
| 5k-10k | `--nproc 2` | 10-15 min | 2-4 GB |
| 10k-50k | `--nproc 1` | 20-60 min | 4-12 GB |
| **54k** | `--nproc 1` | 30-120 min | 8-16 GB |
| >100k | Batch mode needed | hours | Depends on batch size |

---

## What's Changed in This Release

1. **Automatic process limiting** for large datasets
2. **`--nproc` CLI parameter** now functional
3. **Better warnings** about memory and scaling
4. **No code breakage** - backward compatible

---

## Next Steps (Possible Future Work)

- [ ] Implement batch/streaming distance computation
- [ ] Add disk-based caching for intermediate matrices
- [ ] GPU acceleration for distance calculations (RAPIDS)
- [ ] Automatic subsampling recommendation based on available RAM
- [ ] Memory profiling and early-exit if insufficient resources detected

---

## Example Commands

### Test on simulated data (fast validation)
```bash
python alleleatlas.py cgmlst_data/simulated_data.profile test_output --nproc 4
```

### Explore with 5k subset
```bash
python alleleatlas.py cgmlst_5k_subset.csv.gz output_5k --nproc 2
```

### Full run (your dataset, 54k profiles)
```bash
# Option 1: Single process (most stable)
python alleleatlas.py cgmlst_data/cgmlst_573.csv.gz outklebkleb --nproc 1

# Option 2: Dual process (faster, needs more RAM)
python alleleatlas.py cgmlst_data/cgmlst_573.csv.gz outklebkleb --nproc 2 --force
```

### Monitor memory during run (in separate terminal)
```bash
watch -n 2 'ps aux | grep python | grep -v grep | awk "{print \$6}" | paste -sd+ | bc'
```

---

## Questions?

Run with verbose logging:
```bash
python alleleatlas.py cgmlst_data/cgmlst_573.csv.gz outklebkleb --nproc 1 2>&1 | tee run.log
```


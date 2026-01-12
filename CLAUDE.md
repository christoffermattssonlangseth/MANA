# CLAUDE Development Log

This document tracks the development and optimization of the MANA (Metric-Aware Neighborhood Aggregation) methodology, with focus on parameter optimization and aggregation method comparison.

## Project Overview

**MANA** extends CellCharter's neighborhood aggregation with **two-level weighting**:

1. **Hop-level decay**: Neighbors at hop distance `h` are weighted by `hop_decay^h`
2. **Distance-based weighting**: Within each hop, spatial distance determines contribution strength

This creates a continuous, distance-aware spatial context representation that better captures spatial gradients and microenvironmental transitions in spatial transcriptomics data.

## Implementation

### Core Function
`utils/aggregate_neighbors_weighted.py` implements the weighted aggregation approach.

**Key Features:**
- Sparse matrix operations for efficiency
- Multiple distance kernels (exponential, inverse, gaussian, none)
- Flexible hop decay schemes
- Multiple aggregation methods (mean, median, sum, max)
- Sample-aware processing

### Usage Example

```python
from utils.aggregate_neighbors_weighted import aggregate_neighbors_weighted

# Run weighted aggregation
aggregate_neighbors_weighted(
    adata,
    n_layers=3,
    aggregations='mean',
    use_rep='X_scVI',
    out_key='X_weighted_agg_mean',
    hop_decay=0.2,
    distance_kernel='exponential',
    spatial_key='spatial',
    normalize_weights=True,
    include_self=True,
)

# Cluster on weighted features
sc.pp.neighbors(adata, use_rep='X_weighted_agg_mean', n_neighbors=15)
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_weighted')
```

## Parameter Optimization Journey

### Phase 1: Initial Testing (MANA-4)

**Objective:** Systematically test hop_decay and n_layers parameters.

**Parameters Tested:**
- `hop_decay`: [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
- `n_layers`: [2, 3, 4, 5, 6, 7, 8]
- `distance_kernel`: 'exponential' (standard)

**Key Findings:**

| n_layers | Local Purity | Silhouette | Interpretation |
|----------|-------------|------------|----------------|
| 2 | 0.48 | +0.06 | Transcriptionally coherent, spatially weak |
| 3 | 0.54 | +0.04 | Balanced |
| 4 | 0.59 | -0.02 | Transition point |
| 5 | 0.65 | -0.03 | Spatial emphasis |
| 8 | 0.86 | -0.14 | Strong spatial, mixed transcription |

**Best hop_decay:** 0.2 (consistently best across all n_layers)

**Initial Conclusion:** n_layers=8 with hop_decay=0.2 maximized spatial purity (0.857).

### Phase 2: Transcriptional Coherence Analysis (MANA-5 Part 1)

**Critical Discovery:** High spatial purity came at the cost of transcriptional coherence.

**Trade-off Revealed:**
- **High n_layers (7-8)**: Excellent spatial structure BUT negative silhouette scores
  - Clusters are spatially perfect but transcriptionally mixed
  - Not biologically meaningful

- **Low n_layers (2-3)**: Positive silhouette BUT lower spatial purity
  - Transcriptionally coherent
  - Less spatial structure

- **Sweet spot (n_layers=3-4)**: Balanced performance
  - n_layers=3: silhouette=+0.03, purity=0.54
  - n_layers=4: silhouette=-0.02, purity=0.59

**Revised Recommendation:** n_layers=3 (last point with positive silhouette)

### Phase 3: Aggregation Method Comparison (MANA-5 Part 2)

**Objective:** Test different statistical aggregations with optimal parameters.

**Final Parameters Selected:**
- `hop_decay = 0.2` (from MANA-4)
- `n_layers = 3` (from MANA-5 Part 1)
- `distance_kernel = 'exponential'` (standard)

**Aggregation Methods to Test:**
- **mean**: Weighted average (standard)
- **median**: Robust to outliers
- **sum**: Emphasizes total neighborhood activity
- **max**: Highlights strongest signals

**Status:** In progress (MANA-5.ipynb Part 2)

## Evaluation Metrics

### 1. Local Purity (Spatial Coherence)
**Definition:** Fraction of a cell's spatial neighbors that share the same cluster label.

**Formula:**
```
purity(cell_i) = (# neighbors with same label) / (# total neighbors)
local_purity = mean(purity across all cells)
```

**Interpretation:**
- Range: [0, 1]
- Higher = better spatial coherence
- Measures if clusters respect spatial structure

### 2. Silhouette Score (Expression Coherence)
**Definition:** Measures cluster quality in expression space (X_scVI).

**Interpretation:**
- Range: [-1, 1]
- Positive: Cells within clusters are transcriptionally similar
- Negative: Clusters are transcriptionally mixed
- Zero: Overlapping clusters

### 3. Davies-Bouldin Score (Cluster Separation)
**Definition:** Ratio of within-cluster to between-cluster distances.

**Interpretation:**
- Range: [0, ∞)
- Lower = better cluster separation
- Measures compactness and separation in feature space

## Key Insights

### Why Two-Level Weighting?

1. **Hop-level decay addresses network topology:**
   - Not all neighbors are equal
   - Immediate neighbors (hop 1) are more relevant than distant ones (hop 3)
   - Exponential decay: hop_decay^h naturally reflects diminishing influence

2. **Distance weighting addresses spatial heterogeneity:**
   - Cells at same hop distance can be at different spatial distances
   - Network topology ≠ Euclidean distance
   - Closer cells within a hop should contribute more

### Why hop_decay=0.2 Works Best

Lower decay (0.2 vs 0.5) means:
- Hop 1: weight = 0.2
- Hop 2: weight = 0.04
- Hop 3: weight = 0.008

While decay=0.5 gives:
- Hop 1: weight = 0.5
- Hop 2: weight = 0.25
- Hop 3: weight = 0.125

**Lower decay (0.2)** creates steeper falloff → stronger emphasis on immediate neighbors while still allowing distant context to contribute.

### The Spatial-Transcriptional Trade-off

**Fundamental tension:**
- More spatial context (higher n_layers) → better spatial structure but worse transcriptional coherence
- Less spatial context (lower n_layers) → better transcriptional coherence but weaker spatial structure

**Resolution:** n_layers=3-4 balances both objectives.

**Biological interpretation:**
- Positive silhouette = clusters represent transcriptionally distinct cell states
- High local purity = clusters respect tissue spatial organization
- Both are needed for biologically meaningful segmentation

## Workflow Recommendations

### For New Datasets

1. **Build spatial graph:**
   ```python
   import squidpy as sq
   sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
   ```

2. **Test n_layers sweep (decay=0.2, kernel='exponential'):**
   ```python
   for n_layers in [2, 3, 4, 5]:
       aggregate_neighbors_weighted(adata, n_layers=n_layers, hop_decay=0.2)
       # Cluster and evaluate
   ```

3. **Check silhouette scores:**
   - Find last n_layers with positive silhouette
   - This is your optimal depth

4. **Test aggregation methods:**
   - mean (recommended starting point)
   - median (if data is noisy)
   - Compare metrics

5. **Validate biologically:**
   - Marker gene expression
   - Spatial coherence of clusters
   - Alignment with known anatomy

### Quick Start (Recommended Defaults)

```python
# Recommended starting parameters
aggregate_neighbors_weighted(
    adata,
    n_layers=3,
    aggregations='mean',
    use_rep='X_scVI',
    hop_decay=0.2,
    distance_kernel='exponential',
)
```

## File Organization

```
MANA/
├── README.md                    # Project overview
├── CLAUDE.md                    # This file - development log
├── utils/
│   └── aggregate_neighbors_weighted.py  # Core implementation
└── notebooks/
    ├── MANA-4.ipynb            # Initial parameter testing
    └── MANA-5.ipynb            # Optimization and aggregation comparison
```

## Notebooks

### MANA-4.ipynb
**Focus:** Systematic parameter exploration

**Contents:**
- Hop_decay sweep (0.2-0.7)
- n_layers sweep (2-8)
- Local purity analysis
- Biological validation of clusters

**Key Result:** hop_decay=0.2 consistently best

### MANA-5.ipynb
**Focus:** Parameter optimization and aggregation method comparison

**Part 1: n_layers Optimization**
- Trade-off analysis (spatial vs transcriptional coherence)
- Visualization of metric trends
- Recommendation: n_layers=3-4

**Part 2: Aggregation Method Comparison**
- Tests: mean, median, sum, max
- Uses optimal parameters: hop_decay=0.2, n_layers=3
- Comprehensive metric evaluation
- Cluster stability analysis
- Spatial visualization

## Future Directions

### Potential Enhancements

1. **Learnable hop weights:**
   - Instead of fixed hop_decay^h, learn optimal weights per hop
   - Could use cross-validation on local purity

2. **Adaptive distance scales:**
   - Different scales for different tissue regions
   - Adjust based on local cell density

3. **Feature-specific aggregation:**
   - Different aggregation methods for different feature types
   - E.g., mean for some, max for others

4. **Multi-resolution aggregation:**
   - Combine multiple n_layers simultaneously
   - Capture both local and global context

5. **Biological prior integration:**
   - Weight neighbors by cell type compatibility
   - Use ligand-receptor databases

## Technical Notes

### Performance Considerations

- Sparse matrix operations throughout for memory efficiency
- Sample-aware processing to handle batch effects
- Auto-scaling for distance kernels prevents manual tuning

### Edge Cases Handled

- Cells with no neighbors (isolated cells)
- Zero-weight scenarios (normalization)
- Different spatial scales across samples
- Missing spatial coordinates

### Testing Checklist

- [ ] Hop matrices correctly identify cells at exact hop distances
- [ ] Distance weights sum correctly per cell
- [ ] Aggregation methods produce expected outputs
- [ ] Sample-aware processing maintains independence
- [ ] Edge cases don't crash

## References

- **CellCharter:** Original neighborhood aggregation framework
  - Paper: [CellCharter publication]
  - Code: https://github.com/CSOgroup/cellcharter

- **Spatial Transcriptomics:**
  - Squidpy: Spatial single-cell analysis framework
  - Scanpy: Single-cell analysis in Python

## Contact & Contributions

This is an active research project. Parameter recommendations may be refined as we analyze more datasets.

**Current Status:** Phase 3 in progress (aggregation method comparison)

**Last Updated:** 2026-01-12

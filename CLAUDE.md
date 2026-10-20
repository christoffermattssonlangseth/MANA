# CLAUDE Development Log

This document tracks the development and optimization of the MANA (Metric-Aware Neighborhood Aggregation) methodology, with focus on parameter optimization and aggregation method comparison.

## Project Overview

**MANA** extends CellCharter's neighborhood aggregation with **two-level weighting**:

1. **Hop-level decay**: Neighbors at hop distance `h` are weighted by `hop_decay^h`
2. **Distance-based weighting**: Within each hop, spatial distance determines contribution strength

This creates a continuous, distance-aware spatial context representation that better captures spatial gradients and microenvironmental transitions in spatial transcriptomics data.

## 2026-10-20 Re-evaluation (supersedes Phases 1-3 conclusions below)

**Bug:** `hop_decay` was multiplied into the weights *before* row normalization, so it cancelled exactly whenever `normalize_weights=True` (outputs for decay 0.2/0.5/1.0 identical to 1e-16). All MANA-4…MANA-9 results therefore compared only the within-hop distance kernel against CellCharter; "hop_decay=0.2 best" was clustering noise. With `distance_kernel='none'` the function reproduces `cc.gr.aggregate_neighbors` exactly. Fixed: decay now scales each hop's output block; default changed to 1.0. Notebooks that pass `hop_decay=0.2` will now produce different features if re-run.

**Why hop decay is not the lever:** CellCharter concatenates hop rings as separate blocks, so the clusterer learns per-ring importance. A full-covariance GMM (`cc.tl.Cluster`, default) absorbs any fixed block scaling. Decay only matters for scale-sensitive clustering (Leiden on kNN, k-means). The information CellCharter actually discards is *within-ring geometry*.

**Evaluation problems found:**
- `compartment_anno` is a relabelled `cluster_cellcharter_10` (cluster purity 1.000), so it cannot be ground truth for MANA vs CellCharter.
- Local purity / gradient smoothness reward spatial smoothing by construction; the MANA-6 composite score is not evidence of superiority.

**New benchmarks** (scripts were run from a scratch directory; to be moved into the repo if kept):
1. Synthetic tissues with known domains, 7 layouts × 10 seeds, GMM with true K, ARI overall and near boundaries.
2. EAE (877K cells, 107 samples, `..._scVI_embedding_for_mana.h5ad`): ridge regression (sample-grouped CV) from neighbor blocks to anno_L2 composition within 30/60 µm, plus a 1-3 hop composition control.

**Results:**
- Within-hop gaussian weighting (current MANA) helps at holes/tears (boundary ARI 0.36→0.49) and density contrasts (0.30→0.39), is neutral on uniform tissue, and loses on thin layers.
- Variants tested: uniform/kernel blend, per-cell adaptive scale, geodesic distance per hop, pooled multi-scale (geodesic and physical), physical shells, density-adaptive shells, hybrid hops+shells.
- **Winner: physical distance shells** (`aggregate_neighbors_shells`): Gaussian rings at k × median edge length (width 0.5), kept as separate blocks. Synthetic ARI: holes 0.53→0.72, thin layers 0.54→0.61, uniform 0.57→0.65, density contrast 0.63→0.69. EAE 30 µm composition R² 0.376→0.427 (90/107 samples), biggest in white matter (17%), sparse/irregular cells, EAE vs control.
- **Weakness of fixed spacing:** smoothly varying density (ARI 0.65→0.60). Clusters partly follow density (NMI with density bins 0.10 vs 0.02 for CellCharter). Per-cell spacing from noisy Delaunay edges blurred boundaries; hybrid hops+shells doubles dims and hurts GMM.
- **Fix: `adaptive`** — spacing_i = spacing^(1-a) * local^a, local = smoothed mean 15-NN distance. Synthetic mean ARI: a=0 0.630, a=0.5 0.649, a=1 0.644; smoothly varying density 0.60→0.74 (a=0.5). EAE 30 µm R²: a=0 0.427, a=0.5 0.417, a=1 0.395 (full adaptation discards informative density, e.g. sparse white matter; fixed-radius targets also favour a=0). Default a=0.5.

**EAE clustering check** (G3_L1_1, S2_T1_0 control, S2_T1_1, S4_T2_2, S3_B1_0; 63K cells; marker modules from 5101-gene panel):
- K=10: marker eta² equal (0.170 CC / 0.168 shells per cell); same compartments; 17% of cells assigned differently. Marker-based adjudication of disagreeing groups was split (60% of cells favour shells, per-cell 48%).
- 10-seed dorsal-horn recovery: shells 97-99% in 4 samples (CC 91-97%), but shells missed S4_T2_2's dorsal horn in 6/10 seeds (CC 1/10), putting the whole section in its own cluster.
- K=20/30 (5 seeds): CC explains more marker variance (K=30 0.244 vs 0.225, p=0.008); shells slightly more stable (0.72 vs 0.70). 13-15 of 30 clusters are >90% one sample in both.
- `cc.tl.Cluster` is not reproducible even at a fixed random_state on CPU, so method comparisons need several seeds (this also undermines the single-run MANA-6 benchmark).

**S4-T data quality:** S4_T1_1, S4_T1_2, S4_T2_0, S4_T2_2, S4_T3_0 have median 374-685 counts/cell (dataset median 1547) and are the most offset samples in X_scVI (S4_T2_2 rank 1/107; mean rank 12 vs 56). `_scvi_batch` is already per-section (107 batches), so this is lost signal, not missing batch correction. S4_T2_2 is Xenium run 0054723. Neighborhood smoothing makes constant per-section offsets more prominent, so these sections can claim their own clusters.

**Recommended defaults now:** `aggregate_neighbors_shells(adata, n_shells=3, use_rep='X_scVI', sample_key='sample_id')` (adaptive=0.5) then `cc.tl.Cluster`.

## Implementation

### Core Functions

#### 1. Weighted Aggregation
`utils/aggregate_neighbors_weighted.py` implements the weighted aggregation approach.

**Key Features:**
- Sparse matrix operations for efficiency
- Multiple distance kernels (exponential, inverse, gaussian, none)
- Flexible hop decay schemes
- Multiple aggregation methods (mean, median, sum, max, var, std)
- Sample-aware processing

**Usage Example:**

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

#### 2. Spatial Visualization
`utils/plot_spatial_compact_fast.py` provides fast, compact spatial plotting.

**Key Features:**
- Multi-panel layout grouped by sample/condition
- Support for both categorical (clusters) and continuous (gene expression) coloring
- Automatic metadata display (region, course, etc.)
- Category highlighting with custom transparency
- Customizable layouts, colors, and backgrounds
- Publication-ready output with legends/colorbars

**Usage Example:**

```python
from plot_spatial_compact_fast import plot_spatial_compact_fast

# Plot clusters across all samples
plot_spatial_compact_fast(
    adata,
    color='leiden_weighted',      # cluster column or gene name
    groupby='sample_id',           # one panel per sample
    spot_size=8,
    cols=3,                        # 3 columns in grid
    height=8,
    background='white',
    dpi=120
)

# Plot gene expression
plot_spatial_compact_fast(
    adata,
    color='CD8A',                  # gene name
    groupby='sample_id',
    cmap_name='viridis',
    shared_scale=True              # same scale across all panels
)

# Highlight specific clusters
plot_spatial_compact_fast(
    adata,
    color='leiden_weighted',
    groupby='sample_id',
    highlight=['0', '3', '5'],     # highlight clusters 0, 3, 5
    grey_alpha=0.1                 # fade other clusters
)
```

## Parameter Optimization Journey (historical; hop_decay had no effect in these runs)

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

**Aggregation Methods Implemented:**
- **mean**: Weighted average (standard approach, balanced)
- **median**: Weighted median via cumulative weight distribution (robust to outliers)
- **sum**: Weighted sum (emphasizes total neighborhood activity)
- **max**: Maximum weighted contribution (highlights strongest signals)
- **var/std**: Weighted variance/standard deviation (measures neighborhood heterogeneity)

**Implementation Notes:**

*Median Aggregation:*
- Uses weighted quantile approach
- Sorts feature values and computes cumulative weight sum
- Finds 50th percentile using binary search (searchsorted)
- Preserves robustness to outliers while respecting distance weights

*Max Aggregation:*
- Computes weighted contribution = weight × feature value
- Takes maximum across all neighbors per feature
- Useful for detecting dominant microenvironmental signals
- Can be sensitive to noise but highlights strongest influences

**Status:** Ready for testing (MANA-5.ipynb Part 2)

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

### Distance Kernel Choice: Exponential vs Gaussian

**Mathematical Difference:**
- **Exponential**: `weight = exp(-d/scale)` - linear decay in log-space, heavier tail
- **Gaussian**: `weight = exp(-d²/(2*scale²))` - quadratic decay in log-space, sharper cutoff

**Biological Interpretations:**

*Exponential (Current Choice):*
- Models diffusion-like processes (metabolites, cytokines, growth factors)
- Represents gradual microenvironmental influence
- Better for capturing smooth spatial transitions
- Appropriate for: tumor microenvironment gradients, inflammation zones, lesion-associated states

*Gaussian:*
- Models localized, contact-dependent effects
- Represents sharp boundaries between tissue compartments
- Better for discrete anatomical structures
- Appropriate for: organ boundaries, tight cellular niches, well-defined anatomical regions

**Expected Impact on Results:**
- Gaussian would likely produce more clusters with sharper spatial boundaries
- Gaussian would increase local purity but might sacrifice transcriptional coherence
- Exponential better captures the gradual microenvironmental transitions relevant for this project

**Recommendation:** Use exponential kernel as default for microenvironmental analysis. Consider gaussian kernel for datasets with well-defined anatomical compartments.

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

### Quick Start (Recommended Defaults, updated 2026-10-20)

```python
aggregate_neighbors_shells(
    adata,
    n_shells=3,
    use_rep='X_scVI',
    sample_key='sample_id',
)
```

## File Organization

```
MANA/
├── README.md                               # Project overview
├── CLAUDE.md                               # This file - development log
├── utils/
│   ├── aggregate_neighbors_weighted.py    # Weighted aggregation implementation
│   └── plot_spatial_compact_fast.py       # Spatial visualization function
└── notebooks/
    ├── MANA-4.ipynb                       # Initial parameter testing
    ├── MANA-5.ipynb                       # Optimization and aggregation comparison
    └── MANA-6.ipynb                       # Benchmarking vs CellCharter
```

## Notebooks

### MANA-4.ipynb
**Focus:** Systematic parameter exploration

**Contents:**
- Hop_decay sweep (0.2-0.7)
- n_layers sweep (2-8)
- Local purity analysis
- Biological validation of clusters
- Development of `plot_spatial_compact_fast` visualization function

**Key Results:**
- hop_decay=0.2 consistently best
- Custom spatial plotting function created for publication-quality figures

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
- Cluster stability analysis (ARI, NMI)
- Spatial visualization using `plot_spatial_compact_fast`

**Key Features:**
- Module reload mechanism for development
- Imports visualization function from utils
- Ready for aggregation method testing

### MANA-6.ipynb
**Focus:** Head-to-head benchmarking against CellCharter

**Comparison Methods:**
- **CellCharter**: Uniform weights, no distance weighting (hop_decay=1.0, kernel='none')
- **MANA (exponential)**: Distance-weighted with exponential kernel (hop_decay=0.2)
- **MANA (gaussian)**: Distance-weighted with gaussian kernel (hop_decay=0.2)

**Analysis Components:**

1. **Quantitative Metrics**
   - Local purity (spatial coherence)
   - Gradient smoothness (ability to capture transitions)
   - Silhouette score (expression coherence)
   - Davies-Bouldin score (cluster separation)

2. **Statistical Testing**
   - Per-cell purity distributions
   - Wilcoxon signed-rank tests
   - Distribution visualizations (box plots, violin plots)

3. **Cluster Stability**
   - Bootstrap subsampling (10 iterations, 80% samples)
   - ARI scores across iterations
   - Robustness assessment

4. **Biological Validation**
   - Marker gene alignment scores
   - Cell type marker coherence within clusters
   - Customizable marker gene dictionary

5. **Visual Comparison**
   - Side-by-side spatial plots
   - Marker gene expression overlay
   - Multi-sample visualization

6. **Comprehensive Summary**
   - Composite scoring across all metrics
   - Method ranking
   - Use case recommendations

**Key Metrics:**
- `local_purity`: Fraction of neighbors with same cluster label (higher = better)
- `gradient_smoothness`: Average label change across edges (lower = smoother)
- `marker_alignment`: Silhouette on marker genes (higher = better)
- `stability_ari`: Bootstrap consistency (higher = more stable)

**Benchmark Results (MS Lesion Dataset):**
- **MANA (gaussian) wins overall**: Composite score 0.693 vs CellCharter 0.252
- **Local purity**: MANA +14% improvement (0.761 vs 0.665, p < 0.001)
- **Gradient smoothness**: MANA -30% improvement (1.119 vs 1.607, lower = better)
- **Stability**: MANA +5.6% improvement (0.624 vs 0.591 ARI)
- **Silhouette trade-off**: Small decrease (-0.037 vs +0.032) acceptable for spatial gains

**Key Findings:**
- MANA dramatically better at capturing smooth spatial gradients (30% improvement)
- MANA produces significantly more spatially coherent clusters (14% improvement)
- MANA generates more stable/reproducible clusters (5.6% improvement)
- Small transcriptional coherence trade-off consistent with MANA-5 analysis
- Gaussian kernel optimal for tissues with defined boundaries but gradual transitions

**Interpretation Guide:**
- Notebook includes comprehensive results interpretation section
- Ready-to-use methods text for publications
- Biological explanations for why MANA outperforms CellCharter
- Specific recommendations based on tissue characteristics

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

**Current Status (2026-10-20):** Re-evaluated, see "2026-10-20 Re-evaluation" at the top. Distance shells are the recommended method; the benchmark claims below are superseded.

**Recent Updates (2026-01-12):**

*Phase 1 - Parameter Optimization:*
- Added median and max aggregation methods to `aggregate_neighbors_weighted.py`
- Created reusable `plot_spatial_compact_fast.py` visualization function
- Migrated spatial plotting from MANA-4 to utils for reuse across notebooks
- Updated MANA-5 to import and use utils functions
- Added distance kernel comparison theory (exponential vs gaussian)

*Phase 2 - Benchmarking Framework:*
- Created MANA-6.ipynb for comprehensive CellCharter comparison
- Implemented gradient smoothness metric for capturing spatial transitions
- Added marker gene alignment scoring for biological validation
- Built bootstrap stability analysis for cluster robustness testing
- Developed composite scoring system across multiple metrics
- Included statistical significance testing (Wilcoxon signed-rank)
- Added use case recommendations for method selection

*Phase 3 - Benchmark Results & Interpretation:*
- Completed full benchmark on MS lesion dataset (107K cells, 9 samples)
- MANA (gaussian) wins with 0.693 composite score vs 0.252 for CellCharter
- Quantified advantages: +14% spatial coherence, +30% gradient smoothness, +5.6% stability
- Statistical validation: p < 0.001 for spatial improvements
- Added comprehensive interpretation guide to MANA-6 notebook
- Included ready-to-use methods text for publications
- Documented biological rationale for MANA's superiority

**Last Updated:** 2026-01-12

# Metric-Aware Neighborhood Aggregation (MANA)

MANA is a metric-aware generalization of the CellCharter neighborhood aggregation scheme for spatial transcriptomics analysis.

## Overview

[CellCharter](https://github.com/CSOgroup/cellcharter) represents each cell as `[self | mean of 1-hop neighbors | mean of 2-hop | mean of 3-hop]`. The hop rings are kept as separate blocks, so the downstream clusterer already learns how much each ring matters; a fixed hop decay adds nothing on top of that. What CellCharter does lose is geometry: within a ring every cell counts equally, and in irregular tissue one graph hop can span 10 µm or 60 µm.

MANA makes the neighborhood **metric-aware**:

1. **Distance shells** (`aggregate_neighbors_shells`, recommended): rings are defined by physical distance (soft Gaussian rings at ~1, 2, 3 cell spacings) instead of graph hops, still kept as separate blocks.
2. **Within-hop distance weighting** (`aggregate_neighbors_weighted`): keeps graph-hop rings but weights cells within each ring by distance.
3. **Ring variance** (`n_var_shells` / `n_var_layers`): adds the variance of each ring next to its mean, so a mixed neighborhood no longer looks like its average.

Like CellCharter, everything is computed from the cell embedding alone (no annotations).

**What actually matters (2026-10-21 benchmarks):** the biggest single gain is ring variance, and CellCharter already has it: `cc.gr.aggregate_neighbors(..., aggregations=['mean', 'var'])` (off by default). On Allen MERFISH cortex it raises layer ARI from 0.39 to 0.46. MANA shells + variance ties that on regular brain tissue and is better where physical distance matters: irregular synthetic tissue and the cell-type mix within 30-60 µm in EAE. On brain, much of the variance gain is spatial smoothing: smoothing plain CellCharter's labels over the spatial graph (`benchmarks/methods.py: smooth_labels`) recovers most of it and keeps CellCharter's advantage on hippocampal fields. On the synthetic tissues it does not, because there domains differ in their cell-type mix, which only variance captures.

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/MANA.git
cd MANA

# Install dependencies
pip install -r requirements.txt
```

## Quick Start

```python
import scanpy as sc
import squidpy as sq
import cellcharter as cc
from utils import aggregate_neighbors_shells, plot_spatial_compact_fast

# Load your spatial data
adata = sc.read_h5ad('your_data.h5ad')

# Build spatial neighborhood graph (used to estimate the cell spacing)
sq.gr.spatial_neighbors(adata, library_key='sample_id', coord_type='generic', delaunay=True)
cc.gr.remove_long_links(adata)

# Aggregate in physical-distance shells, with the variance of each shell
aggregate_neighbors_shells(
    adata,
    n_shells=3,             # rings at ~1, 2, 3 cell spacings
    n_var_shells=3,         # add per-shell variance blocks
    use_rep='X_scVI',
    sample_key='sample_id',
    out_key='X_mana_shells',
)

# Graph-hop equivalent (= CellCharter with aggregations=['mean', 'var']):
# aggregate_neighbors_weighted(adata, n_layers=3, n_var_layers=3, distance_kernel='none', ...)

# Cluster on aggregated features (as in CellCharter)
model = cc.tl.Cluster(n_clusters=10)
model.fit(adata, use_rep='X_mana_shells')
adata.obs['mana_cluster'] = model.predict(adata, use_rep='X_mana_shells')

# Visualize results
plot_spatial_compact_fast(
    adata,
    color='mana_cluster',
    groupby='sample_id',
    cols=3,
    height=8
)
```

## Repository Structure

```
MANA/
├── README.md                  # This file
├── CLAUDE.md                  # Detailed development log
├── requirements.txt           # Python dependencies
├── benchmarks/                # Synthetic, Allen MERFISH and EAE benchmarks (see benchmarks/README.md)
├── utils/                     # Core implementation
│   ├── __init__.py
│   ├── aggregate_neighbors_weighted.py
│   └── plot_spatial_compact_fast.py
└── notebooks/                 # Analysis notebooks
    ├── README.md
    ├── MANA-4.ipynb          # Parameter optimization
    ├── MANA-5.ipynb          # Aggregation method comparison
    ├── MANA-6.ipynb          # Benchmarking vs CellCharter
    └── archive/              # Older exploratory notebooks
```

## Key Features

### Weighted Aggregation
- Physical-distance shells (recommended)
- Per-ring variance blocks (raw or log)
- Within-hop distance kernels (exponential, inverse, gaussian)
- Optional per-hop block scaling (`hop_decay`, only relevant for Leiden/k-means)
- Various aggregation methods (mean, median, sum, max, var)
- Sparse matrix operations for efficiency

### Spatial Visualization
- Multi-panel layouts grouped by sample/condition
- Categorical and continuous coloring
- Gene expression visualization
- Category highlighting with custom transparency
- Publication-ready output

## Parameters

- **n_shells / n_layers = 3**
- **n_var_shells / n_var_layers = 3**: variance of every ring (raw; `var_log=True` for log variance, which lost less on hippocampal fields but less overall)
- **width = 0.5** (shells), spacing auto-estimated per sample as the median Delaunay edge
- **adaptive = 0.5** (shells): partial adaptation of the spacing to local cell density. Use 0 if physical distance is what matters (e.g. diffusion range); use 1 if density varies independently of the domains you want
- **distance_kernel = 'gaussian'** if using within-hop weighting
- **hop_decay = 1.0** (default). It scales each hop's block. A full-covariance GMM's likelihood is unaffected by this, but `cc.tl.Cluster` starts from k-means, which is scale-sensitive, so block scaling can change results. On realistic data (EAE, realistic synthetic) z-scoring all features made no difference.

## Benchmarking Results (2026-10-21)

All clustering uses `cc.tl.Cluster`, which is not reproducible across runs even at a fixed seed, so every number is a mean over several seeds. Scripts are in [benchmarks/](benchmarks/).

| | Allen cortex layers (K=6, ARI) | Allen layers, whole sections (K=20) | Allen CA1/CA3/DG (K=20) | Realistic synthetic (ARI) | EAE composition within 30 µm (R²) |
|---|---|---|---|---|---|
| CellCharter (default, mean) | 0.393 | 0.408 | **0.253** | 0.224 | 0.554 |
| CellCharter `aggregations=['mean','var']` | 0.457 | **0.489** | 0.129 | 0.297 | 0.577 |
| MANA shells (mean) | 0.393* | 0.397* | 0.237* | 0.260 | 0.609 |
| **MANA shells + variance** | **0.468** | 0.466 | 0.134 | **0.314** | **0.637** |

\* from the first Allen run (cortex 5 seeds, whole sections 3 seeds).

- **Allen MERFISH** (6 coronal sections, 642K cells, labels from CCF registration; cortex 10 seeds, whole sections 5): variance is the gain (layers p<0.001 vs mean only). Shells + variance ties CellCharter mean+var (p=0.15-0.69). Both lose on hippocampal fields, where the mean-only features do best; log variance (`var_log=True`) lost less there (0.29) but gained less on layers (0.44-0.47).
- **Realistic synthetic** (8 layouts x 10 seeds, simulator calibrated to EAE's within-domain structure): shells + variance beats CellCharter mean+var by +0.017 (43/80 runs, p=0.03) and mean-only CellCharter by +0.09; largest at density contrast (0.21 -> 0.49) and small niches.
- **Label smoothing and consensus** (annotation-free post-processing: 2 rounds of majority vote over each cell's Delaunay neighbours; Hungarian-aligned majority vote over 5 seeds). Allen K=20, layers / CA1-CA3-DG / regions: CellCharter 0.408 / 0.253 / 0.305 -> smoothed 0.491 / 0.259 / 0.339 -> consensus + smoothed **0.513 / 0.305 / 0.364**, the best on all three; CellCharter mean+var 0.489 / 0.129 -> smoothed 0.508 / 0.132; shells + variance 0.466 / 0.134 -> 0.472 / 0.139. Smoothing lifts CellCharter's layer ARI in every seed (+0.07 to +0.11) but variance features only by 0-0.02, as their labels are already smooth. Realistic synthetic, smoothed: CellCharter 0.259, mean+var 0.311, shells 0.278, shells + variance 0.323 (vs mean+var +0.012, p=0.16; unsmoothed in the same run +0.021, p=0.007).
- **Higher K does not resolve the layers vs hippocampus trade-off.** At K=35 the variance features still score 0.14-0.16 on CA1/CA3/DG vs 0.24-0.26 for CellCharter; only thalamic nuclei improve (0.10 -> 0.16).
- **Negative results:** neighbour cell-type composition from hard k-means types (40 types) lost on layers (0.359 cortex-only) though it gained on hippocampus (0.315); PCA-50 compression of mean+var features hurt layers (0.296; 0.357 after z-scoring) while recovering hippocampus and thalamus. Hard type bins are also a step away from the embedding-only design.
- **EAE composition** (877K cells, 107 sections; ridge regression from neighbor blocks to the anno_L2 mix around each cell, sample-held-out CV): shells + variance is best for physical radii (30 µm 0.637, 60 µm 0.768, all 107 samples better than CellCharter). CellCharter mean+var is best on the 1-3 hop composition control (0.790), as expected for a graph-hop target.
- **EAE clustering** (5 sections, marker modules, 5 seeds): at K=10 adding variance helps slightly (marker eta² 0.170 -> 0.176-0.179) and shells + variance is the most seed-stable (0.77 vs 0.67); at K=20 all methods tie. The S4-T sections (2-4x lower counts per cell) can claim their own clusters with any method; consider fitting without them and predicting them afterwards.

**Corrections to earlier results:**
- `hop_decay` had no effect before the 2026-10-20 fix, and local purity / gradient smoothness reward smoothing by construction, so the MANA-4/MANA-6 conclusions (`hop_decay=0.2` best; composite 0.693 vs 0.252) do not hold. `compartment_anno` is a relabelled CellCharter (k=10) clustering and cannot be ground truth.
- The first EAE composition numbers (CellCharter 0.376, shells 0.427 at 30 µm) were too low: per-sample subsetting dropped unused cell-type categories and misaligned the targets in 26/107 sections (fixed in `real_composition.py`). Method rankings were unchanged; the adaptive-spacing comparison (a=0 / 0.5 / 1) has not been re-run with the fix.
- The original synthetic simulator draws cells independently within domains, which exaggerates effects (z-scoring features looked like +0.11 ARI there, ~0 on EAE). Use `realistic_synthetic.py`.

## Documentation

- **[CLAUDE.md](CLAUDE.md)**: Comprehensive development log with parameter optimization insights
- **[notebooks/README.md](notebooks/README.md)**: Notebook descriptions and usage
- **[utils/](utils/)**: Function docstrings with detailed parameter descriptions

## Citation

If you use MANA in your research, please cite:

```
[Citation information to be added]
```

## References

- **CellCharter**: Original neighborhood aggregation framework ([GitHub](https://github.com/CSOgroup/cellcharter))
- **Squidpy**: Spatial single-cell analysis ([Documentation](https://squidpy.readthedocs.io/))
- **Scanpy**: Single-cell analysis in Python ([Documentation](https://scanpy.readthedocs.io/))

## Status

**In development.** Ring variance (already available in CellCharter) is the largest feature-side gain; MANA shells + variance matches it on regular brain tissue and improves on it where physical distance matters (irregular tissue, physical-radius composition). Whichever features you use, cluster with several seeds and consider consensus plus neighbour label smoothing: on Allen brain, plain CellCharter with both is the best overall. Next: understand why variance features lose hippocampal fields, and test smoothing/consensus on EAE.
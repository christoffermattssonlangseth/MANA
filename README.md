# Metric-Aware Neighborhood Aggregation (MANA)

MANA is a metric-aware generalization of the CellCharter neighborhood aggregation scheme for spatial transcriptomics analysis.

## Overview

[CellCharter](https://github.com/CSOgroup/cellcharter) represents each cell as `[self | mean of 1-hop neighbors | mean of 2-hop | mean of 3-hop]`. The hop rings are kept as separate blocks, so the downstream clusterer already learns how much each ring matters; a fixed hop decay adds nothing on top of that. What CellCharter does lose is geometry: within a ring every cell counts equally, and in irregular tissue one graph hop can span 10 µm or 60 µm.

MANA makes the neighborhood **metric-aware**:

1. **Distance shells** (`aggregate_neighbors_shells`, recommended): rings are defined by physical distance (soft Gaussian rings at ~1, 2, 3 cell spacings) instead of graph hops, still kept as separate blocks.
2. **Within-hop distance weighting** (`aggregate_neighbors_weighted`): keeps graph-hop rings but weights cells within each ring by distance.

The gains concentrate where tissue geometry is irregular: tears and holes, density changes, white matter and disrupted (EAE) tissue.

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

# Aggregate in physical-distance shells
aggregate_neighbors_shells(
    adata,
    n_shells=3,             # rings at ~1, 2, 3 cell spacings
    use_rep='X_scVI',
    sample_key='sample_id',
    out_key='X_mana_shells',
)

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
- **width = 0.5** (shells), spacing auto-estimated per sample as the median Delaunay edge
- **adaptive = 0.5** (shells): partial adaptation of the spacing to local cell density. Use 0 if physical distance is what matters (e.g. diffusion range); use 1 if density varies independently of the domains you want
- **distance_kernel = 'gaussian'** if using within-hop weighting
- **hop_decay = 1.0** (default). It scales each hop's block; a full-covariance GMM absorbs this, so only tune it for Leiden/k-means pipelines.

## Benchmarking Results (2026-10-20)

**Synthetic tissues with known domains** (10 seeds each, GMM with true K, overall ARI):

| Tissue layout | CellCharter | MANA within-hop (gaussian) | Shells, adaptive=0 | Shells, adaptive=0.5 |
|---|---|---|---|---|
| Holes / tears | 0.53 | 0.57 | 0.72 | 0.72 |
| Thin layers | 0.54 | 0.52 | 0.61 | 0.61 |
| Uniform | 0.57 | 0.56 | **0.65** | 0.63 |
| Density contrast | 0.63 | 0.64 | 0.69 | 0.70 |
| Small niches | 0.50 | 0.50 | 0.51 | 0.51 |
| Smoothly varying density | 0.65 | 0.60 | 0.60 | **0.73** |
| Mean | 0.57 | 0.57 | 0.63 | **0.65** |

The default (adaptive=0.5) beats CellCharter in all six layouts (wins 7-9/10 seeds each; p ≤ 0.014 for holes, thin layers and smoothly varying density).

**EAE spinal cord** (877K cells, 107 samples): how well do the neighbor blocks predict the true cell-type composition (anno_L2) around each cell? (Ridge regression, sample-held-out CV.)

| R² | CellCharter | MANA within-hop | Shells, adaptive=0 | Shells, adaptive=0.5 |
|---|---|---|---|---|
| composition within 30 µm | 0.376 | 0.406 | **0.427** | 0.417 |
| composition within 60 µm | 0.454 | 0.462 | **0.478** | 0.477 |
| 1-3 hop composition (CellCharter-favoured control) | **0.493** | 0.487 | 0.465 | 0.474 |

Shells reduce error most in white matter (12-17%), the sparsest/most irregular cells and EAE vs control tissue (7-9% vs ~0%). The fixed-radius targets favour adaptive=0 by construction.

**EAE clustering** (5 samples, 63K cells, `cc.tl.Cluster`, multiple seeds): at K=10 both methods find the same compartments and explain the same marker-module variance (0.17 per cell). At K=20/30 CellCharter explains slightly more (K=30: 0.244 vs 0.225, all 5 seeds) while MANA shells is slightly more seed-stable (ARI 0.72 vs 0.70). So MANA's better neighborhood encoding does not translate into better clusters on this tissue; its advantage is in irregular geometry and in using the features directly. Single clustering runs are not reproducible enough to compare methods: always use several seeds.

**Data note:** the S4-T sections (S4_T1_1, S4_T1_2, S4_T2_0, S4_T2_2, S4_T3_0) have 2-4x lower counts per cell (median 374-685 vs 1547) and are the most offset samples in the scVI embedding, despite per-section batch correction. They can claim their own clusters with either method; consider fitting clusters without them and predicting them afterwards.

**Known limitation:** with adaptive=0, where density varies smoothly and independently of domains, shells partly cluster by density; adaptive>0 fixes this at a small cost in physical-distance fidelity.

**Superseded results:** the earlier MANA-4/MANA-6 conclusions (`hop_decay=0.2` best; composite score 0.693 vs 0.252) do not hold. `hop_decay` had no effect before the 2026-10-20 fix, and local purity / gradient smoothness reward spatial smoothing by construction. `compartment_anno` is a relabelled CellCharter (k=10) clustering and cannot serve as ground truth for MANA vs CellCharter.

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

**In development.** Distance shells beat CellCharter on geometrically irregular tissue (synthetic ground truth and EAE composition encoding). On EAE domain clustering the two methods are roughly equivalent (CellCharter slightly ahead at high K). Next: a public dataset with expert domain annotations.
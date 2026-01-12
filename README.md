# Metric-Aware Neighborhood Aggregation (MANA)

MANA is a metric-aware generalization of the CellCharter neighborhood aggregation scheme for spatial transcriptomics analysis.

## Overview

The central idea is to extend [CellCharter](https://github.com/CSOgroup/cellcharter) with a **two-level weighting strategy**:

1. **Hop-level decay**: Neighbors at hop distance `h` are weighted by `hop_decay^h`
2. **Distance-based weighting**: Within each hop, spatial distance determines contribution strength

Instead of treating all neighboring cells equally within discrete hop-based rings, MANA weights each neighbor by its physical distance, such that closer cells contribute more strongly than distant ones.

This creates a **continuous, distance-aware spatial context representation** that better captures:
- Spatial gradients and microenvironmental transitions
- Lesion-associated tissue states
- Smooth variations across tissue space

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
from utils import aggregate_neighbors_weighted, plot_spatial_compact_fast

# Load your spatial data
adata = sc.read_h5ad('your_data.h5ad')

# Build spatial neighborhood graph
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)

# Run weighted aggregation
aggregate_neighbors_weighted(
    adata,
    n_layers=3,           # Number of neighborhood layers
    hop_decay=0.2,        # Weight decay per hop
    aggregations='mean',  # Aggregation method
    use_rep='X_scVI',     # Input representation
    distance_kernel='exponential'
)

# Cluster on aggregated features
sc.pp.neighbors(adata, use_rep='X_weighted_agg_mean', n_neighbors=15)
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_weighted')

# Visualize results
plot_spatial_compact_fast(
    adata,
    color='leiden_weighted',
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
    └── archive/              # Older exploratory notebooks
```

## Key Features

### Weighted Aggregation
- Multiple distance kernels (exponential, inverse, gaussian)
- Flexible hop decay schemes
- Various aggregation methods (mean, median, sum, max, var)
- Sparse matrix operations for efficiency

### Spatial Visualization
- Multi-panel layouts grouped by sample/condition
- Categorical and continuous coloring
- Gene expression visualization
- Category highlighting with custom transparency
- Publication-ready output

## Optimal Parameters

Based on systematic testing in [notebooks/MANA-4.ipynb](notebooks/MANA-4.ipynb):

- **hop_decay = 0.2**: Consistently best across all configurations
- **n_layers = 3-4**: Optimal balance between spatial purity and transcriptional coherence
- **distance_kernel = 'exponential'**: Best for gradual microenvironmental transitions

See [CLAUDE.md](CLAUDE.md) for detailed parameter optimization analysis.

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

Active development. Parameter recommendations based on systematic testing are stable and ready for use.
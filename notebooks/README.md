# MANA Notebooks

This directory contains analysis notebooks for the MANA (Metric-Aware Neighborhood Aggregation) project.

## Active Notebooks

### [MANA-4.ipynb](MANA-4.ipynb)
**Parameter Optimization**

Systematic testing of hop_decay and n_layers parameters to find optimal balance between spatial purity and transcriptional coherence.

**Key Findings:**
- `hop_decay = 0.2` performs best across all configurations
- `n_layers = 3-4` provides optimal balance
- Higher n_layers improve spatial structure but reduce transcriptional coherence

### [MANA-5.ipynb](MANA-5.ipynb)
**Aggregation Method Comparison**

Compares different statistical aggregation methods (mean, median, sum, max) using optimal parameters from MANA-4.

**Methods Tested:**
- Mean: Standard weighted average
- Median: Robust to outliers
- Sum: Emphasizes total neighborhood activity
- Max: Highlights strongest signals

## Usage

All notebooks import utility functions from the `utils` directory:

```python
import sys
sys.path.insert(0, '../utils')
from plot_spatial_compact_fast import plot_spatial_compact_fast
from aggregate_neighbors_weighted import aggregate_neighbors_weighted
```

Or alternatively:

```python
import sys
sys.path.insert(0, '..')
from utils import plot_spatial_compact_fast, aggregate_neighbors_weighted
```

## Archive

Older exploratory notebooks are available in the [`archive/`](archive/) directory.

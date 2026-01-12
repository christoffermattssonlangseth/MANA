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

### [MANA-6.ipynb](MANA-6.ipynb)
**Benchmarking Against CellCharter**

Comprehensive head-to-head comparison between MANA and base CellCharter method.

**Comparison Focus:**
- Quantitative metrics (spatial coherence, expression coherence, cluster separation)
- Statistical significance testing (Wilcoxon, distribution analysis)
- Cluster stability via bootstrap subsampling
- Biological validation with marker gene alignment
- Visual side-by-side comparison
- Gradient smoothness analysis

**Methods Compared:**
- CellCharter (uniform weights, no distance weighting)
- MANA with exponential kernel (hop_decay=0.2)
- MANA with gaussian kernel (hop_decay=0.2)

**Key Metrics:**
- `local_purity`: Spatial coherence
- `gradient_smoothness`: Ability to capture smooth transitions
- `marker_alignment`: Biological validation
- `stability_ari`: Robustness across subsampling

**Benchmark Results (MS Lesion Dataset, 107K cells):**
- 🏆 MANA (gaussian) wins overall: 0.693 vs 0.252 composite score
- ✅ +14% spatial coherence (p < 0.001)
- ✅ +30% better gradient smoothness
- ✅ +5.6% more stable clusters
- ⚠️ Small transcriptional trade-off acceptable

**Interpretation Guide:**
The notebook includes a comprehensive "Results Interpretation Guide" section with:
- Detailed metric-by-metric analysis with winners
- Biological interpretations (why MANA works better)
- Key claims for publications with exact percentages
- Ready-to-use methods text for papers
- Recommendations based on tissue characteristics
- Composite score breakdown and explanation

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

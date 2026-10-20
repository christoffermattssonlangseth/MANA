# MANA Benchmarks (2026-10-20)

Scripts behind the benchmark tables in [../README.md](../README.md). Run them with the `cellcharter` conda env.

## Synthetic tissues with known domains

`synthetic_benchmark.py` simulates seven tissue layouts (uniform, density_contrast, thin_layers,
small_niches, holes, variable_density, gradient). Domains are defined by cell-type composition and
cells get noisy latent profiles. Methods are clustered with a full-covariance GMM (and k-means) at
the true K and scored with ARI overall, near boundaries (≤2 hops from a different domain) and in
the interior.

```bash
python synthetic_benchmark.py 10 2,3 results.csv     # seeds, n_layers, output
```

`variants.py` holds the candidate aggregation schemes (blend, adaptive, geo_hop, geo/phys
multiscale, phys_shells, alpha<a>, hybrid, repo_shells_<a>). `run_variants.py` benchmarks any of
them against CellCharter and current MANA (GMM, n_layers=3, 10 seeds):

```bash
python run_variants.py repo_shells_0,repo_shells_0.5 results_variants.csv
```

## EAE composition encoding

`real_composition.py` asks how well each representation's neighbor blocks (self block excluded)
predict the anno_L2 composition around each cell within 30 / 60 µm, plus a 1-3 hop composition
control that favours CellCharter. Ridge regression with GroupKFold over samples. Gains are
stratified by local geometry, density, compartment, cell class and condition.

```bash
python real_composition.py 3 repo_shells_0,repo_shells_0.5 eae_mana_subset.h5ad
```

The input is a light extract (obs annotations + `spatial` + `X_scVI`, ~130 MB) of
`/Volumes/processing2/RRmap/data/EAE_proseg_clustered_louvain_leiden_all_sections_annotated_rotated_scVI_embedding_for_mana.h5ad`,
read with h5py so the full file is never loaded:

```python
import h5py, anndata as ad, pandas as pd, numpy as np
from anndata.io import read_elem
f = h5py.File(PATH, 'r')
cols = ['sample_id', 'sample_name', 'region', 'condition', 'course',
        'anno_L1', 'anno_L2', 'anno_L3', 'compartment_anno', 'volume']
obs = pd.DataFrame({c: read_elem(f['obs'][c]) for c in cols},
                   index=read_elem(f['obs'][f['obs'].attrs['_index']]))
a = ad.AnnData(obs=obs, obsm={'spatial': f['obsm']['spatial'][:],
                              'X_scVI': f['obsm']['X_scVI'][:].astype(np.float32)})
a.write_h5ad('eae_mana_subset.h5ad', compression='gzip')
```

Note: `compartment_anno` is a relabelled `cluster_cellcharter_10` and is only used for stratifying
results, never as ground truth.

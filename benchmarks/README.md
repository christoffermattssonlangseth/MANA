# MANA Benchmarks (2026-10-21)

Scripts behind the benchmark tables in [../README.md](../README.md). Run them with the `cellcharter` conda env.

## Harness

`methods.py` is the registry of feature builders (CellCharter, CellCharter mean+var, MANA shells,
MANA shells+var) plus shared graph building and `cc.tl.Cluster` helpers. Add a new idea to
`METHODS` and the realistic synthetic and Allen benchmarks pick it up. All features must be
computed from the embedding alone (no annotations); annotations are only used for scoring.
`cc.tl.Cluster` is not reproducible across runs, so every benchmark uses several seeds.
`smooth_labels` (majority vote over graph neighbours) and `consensus_labels` (Hungarian-aligned
vote over seeds) are annotation-free post-processing steps; `realistic_synthetic.py` reports ARI
both raw and smoothed.

## Realistic synthetic tissues (use this one)

`realistic_synthetic.py` adds smooth within-domain cell-state variation and cell-type patches,
calibrated to EAE's block-variance profile, plus lesion-rim and patchy-vs-mixed layouts.
8 layouts x N seeds, `cc.tl.Cluster` at the true K, ARI.

```bash
N_JOBS=8 python realistic_synthetic.py 10 results_realistic.csv
```

## Allen MERFISH whole mouse brain (independent atlas labels)

`allen_merfish.py` downloads 6 coronal sections of MERFISH-C57BL6J-638850 (~1.6 GB, public S3),
builds PCA + per-section features for every method in `methods.py`, and scores `cc.tl.Cluster`
against CCF labels: cortical layers (cortex-only K=6, 10 seeds) and layers / CA1-CA3-DG /
thalamic nuclei / regions in whole-section clusterings (K=20, 5 seeds). Clustering 642K cells is
heavy; run it on a big machine.

```bash
python allen_merfish.py download /path/allen && python allen_merfish.py prep /path/allen
N_JOBS=5 python allen_merfish.py bench /path/allen
```

## Synthetic tissues with known domains (original, i.i.d. cells)

Cells are drawn independently within domains, which exaggerates some effects (e.g. feature
z-scoring); prefer `realistic_synthetic.py`.

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

Fixed 2026-10-21: one-hot targets were built from per-sample category codes, and AnnData
subsetting drops unused categories, so targets were misaligned in 26/107 sections. Results from
before the fix have too-low R² (rankings were unchanged).

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

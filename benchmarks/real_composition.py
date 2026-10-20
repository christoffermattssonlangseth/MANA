"""Real-data test on EAE: how well does each aggregated neighbourhood representation encode
the true local cell-type composition (anno_L2, expression-derived, independent of both methods)?

Targets (per cell, excluding self):
  phys30 / phys60 : anno_L2 composition among cells within 30 / 60 um      (physical niche)
  hop3            : anno_L2 composition among 1..3-hop graph neighbours     (CellCharter-favoured control)
Features: neighbour blocks only (hops 1..L), so the cell's own profile is not used.
Model: ridge, GroupKFold over samples.
"""
import os, sys, io, contextlib, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np, pandas as pd, anndata as ad, squidpy as sq, cellcharter as cc
from scipy.sparse import csr_matrix, eye
from scipy.spatial import cKDTree
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from utils.aggregate_neighbors_weighted import aggregate_neighbors_weighted
from variants import features as vfeatures
KINDS = sys.argv[2].split(',') if len(sys.argv) > 2 else ['geo_hop']

L = int(sys.argv[1]) if len(sys.argv) > 1 else 3
a = ad.read_h5ad(sys.argv[3] if len(sys.argv) > 3 else 'eae_mana_subset.h5ad')
types = a.obs.anno_L2.astype('category').cat.categories
LAT = a.obsm['X_scVI'].shape[1]

feats = {m: [] for m in ['CellCharter', 'MANA-gauss'] + KINDS}
targets = {t: [] for t in ['phys30', 'phys60', 'hop3']}
meta = []
for s, idx in a.obs.groupby('sample_id', observed=True).indices.items():
    b = a[idx].copy()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        sq.gr.spatial_neighbors(b, coord_type='generic', delaunay=True); cc.gr.remove_long_links(b)
        cc.gr.aggregate_neighbors(b, n_layers=L, use_rep='X_scVI', out_key='cc')
        for k in ['gaussian']:
            aggregate_neighbors_weighted(b, n_layers=L, use_rep='X_scVI', out_key=k, hop_decay=0.2,
                                         distance_kernel=k, normalize_weights=True)
    feats['CellCharter'].append(b.obsm['cc'][:, LAT:])
    for k in KINDS:
        feats[k].append(vfeatures(b, L, k, use_rep='X_scVI')[:, LAT:])
    feats['MANA-gauss'].append(b.obsm['gaussian'][:, LAT:])

    onehot = np.eye(len(types))[b.obs.anno_L2.cat.codes.values]
    xy = b.obsm['spatial']; tree = cKDTree(xy)
    for r, name in [(30, 'phys30'), (60, 'phys60')]:
        P = tree.sparse_distance_matrix(tree, r, output_type='coo_matrix').tocsr()
        P.data[:] = 1; P = P - eye(len(xy), format='csr').multiply(P); P.eliminate_zeros()
        n = np.asarray(P.sum(1)).ravel(); targets[name].append((P @ onehot) / np.maximum(n, 1)[:, None])
    A = b.obsp['spatial_connectivities'].tocsr(); A.data[:] = 1
    R = A.copy(); F = A.copy()
    for _ in range(L - 1):
        F = F @ A; F.data[:] = 1; R = R + F; R.data[:] = 1
    R = R - eye(len(xy), format='csr').multiply(R); R.eliminate_zeros()
    n = np.asarray(R.sum(1)).ravel(); targets['hop3'].append((R @ onehot) / np.maximum(n, 1)[:, None])

    D = b.obsp['spatial_distances'].tocsr()
    mx = D.max(1).toarray().ravel(); mn = np.array([D.data[D.indptr[i]:D.indptr[i + 1]].min() if D.indptr[i + 1] > D.indptr[i] else np.nan for i in range(D.shape[0])])
    dens = np.array([len(x) for x in tree.query_ball_point(xy, 50)])
    meta.append(pd.DataFrame(dict(sample=s, ratio=mx / mn, density50=dens, deg=np.diff(A.indptr),
                                  comp=b.obs.compartment_anno.values, L1=b.obs.anno_L1.values,
                                  condition=b.obs.condition.values, course=b.obs.course.values), index=b.obs_names))

meta = pd.concat(meta)
groups = meta['sample'].astype('category').cat.codes.values
feats = {k: np.vstack(v) for k, v in feats.items()}
targets = {k: np.vstack(v) for k, v in targets.items()}
gkf = GroupKFold(5)
summary, percell = [], {}
for t, Y in targets.items():
    for m, X in feats.items():
        pred = np.zeros_like(Y)
        for tr, te in gkf.split(X, Y, groups):
            mu, sd = X[tr].mean(0), X[tr].std(0) + 1e-9
            pred[te] = Ridge(alpha=1.0).fit((X[tr] - mu) / sd, Y[tr]).predict((X[te] - mu) / sd)
        sse = ((pred - Y) ** 2).sum(1)
        percell[(t, m)] = sse
        r2 = 1 - ((pred - Y) ** 2).sum() / ((Y - Y.mean(0)) ** 2).sum()
        per_sample = pd.Series(sse).groupby(groups).mean()
        summary.append(dict(target=t, method=m, R2=r2, per_sample_sse=per_sample))
        print(f'{t:7s} {m:12s} R2={r2:.4f}', flush=True)

# paired per-sample comparison vs CellCharter
from scipy.stats import wilcoxon
rows = []
for t in targets:
    base = [s for s in summary if s['target'] == t and s['method'] == 'CellCharter'][0]
    for s in summary:
        if s['target'] != t or s['method'] == 'CellCharter':
            continue
        d = base['per_sample_sse'] - s['per_sample_sse']  # >0 means MANA better
        rows.append(dict(target=t, method=s['method'], R2_CC=base['R2'], R2_MANA=s['R2'],
                         samples_won=f'{(d > 0).sum()}/{len(d)}', p=wilcoxon(d).pvalue))
print(pd.DataFrame(rows).round(4).to_string(index=False))

# where does MANA-gauss help? per-cell error reduction on phys30, stratified
for t in ['phys30', 'phys60']:
    for m in ['MANA-gauss'] + KINDS:
        meta[f'{m}|{t}'] = (percell[(t, 'CellCharter')] - percell[(t, m)]) / percell[(t, 'CellCharter')].mean()
meta['ratio_bin'] = pd.qcut(meta['ratio'], 5).astype(str)
meta['density_bin'] = pd.qcut(meta['density50'], 5, duplicates='drop').astype(str)
for col in ['ratio_bin', 'density_bin', 'comp', 'L1', 'condition']:
    print(f'\n--- relative error reduction by {col} (positive = MANA-gauss better)')
    print(meta.groupby(col, observed=True)[[c for c in meta.columns if '|' in c]].mean().round(4).assign(
        n=meta.groupby(col, observed=True).size()).to_string())
meta.to_parquet(f'real2_percell_L{L}.parquet')
pd.DataFrame(rows).to_csv(f'real2_summary_L{L}.csv', index=False)

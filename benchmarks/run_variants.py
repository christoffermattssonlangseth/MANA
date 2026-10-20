import sys, io, contextlib, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, squidpy as sq, cellcharter as cc
from sklearn.mixture import GaussianMixture
from sklearn.metrics import adjusted_rand_score
import synthetic_benchmark as bench
from variants import features
from utils.aggregate_neighbors_weighted import aggregate_neighbors_weighted

KINDS = sys.argv[1].split(',') if len(sys.argv) > 1 else ['repo_shells_0.5']


def run_one(name, seed, L=3):
    a, cont = bench.simulate(name, seed)
    sq.gr.spatial_neighbors(a, coord_type='generic', delaunay=True); cc.gr.remove_long_links(a)
    truth = a.obs['truth'].values; K = len(np.unique(truth)); bm = bench.boundary_mask(a)
    cc.gr.aggregate_neighbors(a, n_layers=L, use_rep='X_lat', out_key='cc')
    aggregate_neighbors_weighted(a, n_layers=L, use_rep='X_lat', out_key='g', distance_kernel='gaussian')
    F = {'CellCharter': a.obsm['cc'], 'MANA-gauss': a.obsm['g']}
    for k in KINDS:
        F[k] = features(a, L, k)
    rows = []
    for m, X in F.items():
        pred = GaussianMixture(K, covariance_type='full', n_init=3, random_state=seed).fit_predict(X)
        row = dict(scenario=name, seed=seed, method=m, ARI=adjusted_rand_score(truth, pred),
                   ARI_boundary=adjusted_rand_score(truth[bm], pred[bm]))
        if name == 'small_niches': row['niche_F1'] = bench.matched_f1(truth, pred, 2)
        if name == 'thin_layers': row['thin_F1'] = bench.matched_f1(truth, pred, 1)
        rows.append(row)
    return rows


if __name__ == '__main__':
    out = []
    for name in [s for s in bench.SCENARIOS if s != 'gradient']:
        for s in range(10):
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                out += run_one(name, s)
        print('done', name, flush=True)
        pd.DataFrame(out).to_csv(sys.argv[2] if len(sys.argv) > 2 else 'results_variants.csv', index=False)

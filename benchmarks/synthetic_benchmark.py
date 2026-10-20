"""Ground-truth synthetic benchmark: CellCharter vs MANA distance-weighted aggregation.

Each scenario simulates a tissue with known spatial domains. Domains are defined by
cell-type *composition* (like real niches), cells get a noisy latent profile from their
type. We aggregate with each method, cluster with the true K, and score against truth.
"""
import os, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np, pandas as pd, anndata as ad, squidpy as sq, cellcharter as cc
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, f1_score
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import cross_val_predict
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree
from utils.aggregate_neighbors_weighted import aggregate_neighbors_weighted

N_TYPES, LAT = 8, 10


def type_profiles(rng):
    return rng.normal(0, 1.0, (N_TYPES, LAT))


def compositions(rng, k, conc=0.6):
    # distinct but overlapping niche compositions
    return rng.dirichlet(np.full(N_TYPES, conc), size=k)


def poisson(rng, n, box=(0, 100, 0, 100)):
    x0, x1, y0, y1 = box
    return np.c_[rng.uniform(x0, x1, n), rng.uniform(y0, y1, n)]


def jitter_grid(rng, spacing, box=(0, 100, 0, 100), jit=0.35):
    xs = np.arange(box[0], box[1], spacing); ys = np.arange(box[2], box[3], spacing)
    g = np.array(np.meshgrid(xs, ys)).reshape(2, -1).T
    return g + rng.normal(0, jit * spacing, g.shape)


# ---------------- scenarios: return coords, domain labels (int), optional continuous truth
def sc_uniform(rng):
    xy = jitter_grid(rng, 1.3)
    seeds = rng.uniform(0, 100, (6, 2))
    lab = cKDTree(seeds).query(xy)[1]
    return xy, lab, None


def sc_density_contrast(rng):
    # left half dense (like hypercellular lesion rim), right half sparse; domains are vertical bands
    dense = jitter_grid(rng, 0.9, (0, 50, 0, 100)); sparse = jitter_grid(rng, 2.2, (50, 100, 0, 100))
    xy = np.vstack([dense, sparse])
    lab = np.digitize(xy[:, 0], [25, 50, 75])
    return xy, lab, None


def sc_thin_layers(rng):
    # cortex-like stacked layers of very different thickness
    xy = jitter_grid(rng, 1.2)
    edges = np.cumsum([30, 4, 20, 3, 25])  # thin layers of 4 and 3 units (~3 cells)
    lab = np.digitize(xy[:, 1], edges)
    lab = np.array([0, 1, 2, 1, 3, 4])[lab]  # thin layers share a domain label
    return xy, lab, None


def sc_small_niches(rng):
    # background + many small perivascular-cuff-like niches (radius 2.5-4)
    xy = jitter_grid(rng, 1.2)
    centers = rng.uniform(8, 92, (25, 2)); radii = rng.uniform(2.5, 4.0, 25)
    d = np.linalg.norm(xy[:, None] - centers[None], axis=2)
    lab = (d < radii).any(1).astype(int)
    # a second large background split, so K=3
    lab = np.where(lab == 1, 2, (xy[:, 0] > 50).astype(int))
    return xy, lab, None


def sc_holes(rng):
    # tissue with tears/holes: Delaunay links span gaps (partially removed by remove_long_links)
    xy = jitter_grid(rng, 1.2)
    keep = np.ones(len(xy), bool)
    for _ in range(12):
        c = rng.uniform(5, 95, 2); r = rng.uniform(3, 7)
        keep &= np.linalg.norm(xy - c, axis=1) > r
    # long narrow cracks along domain boundaries
    keep &= ~((np.abs(xy[:, 0] - 50) < 1.6) & (xy[:, 1] > 15))
    keep &= ~((np.abs(xy[:, 1] - 50) < 1.6) & (xy[:, 0] < 85))
    xy = xy[keep]
    lab = (xy[:, 0] > 50).astype(int) * 2 + (xy[:, 1] > 50).astype(int)
    return xy, lab, None


def sc_variable_density(rng):
    # density varies smoothly AND independently of domain (realistic heterogeneity)
    n = 9000
    xy = poisson(rng, n * 3)
    p = 0.15 + 0.85 * (0.5 + 0.5 * np.sin(xy[:, 0] / 9) * np.cos(xy[:, 1] / 13))
    xy = xy[rng.uniform(size=len(xy)) < p][:n]
    seeds = rng.uniform(0, 100, (6, 2))
    lab = cKDTree(seeds).query(xy)[1]
    return xy, lab, None


def sc_gradient(rng):
    # continuous composition gradient (e.g. lesion core -> periplaque); truth = continuous
    xy = jitter_grid(rng, 1.2)
    r = np.linalg.norm(xy - 50, axis=1)
    t = np.clip(r / 50, 0, 1)
    lab = np.digitize(t, [0.25, 0.5, 0.75])
    return xy, lab, t


SCENARIOS = dict(uniform=sc_uniform, density_contrast=sc_density_contrast, thin_layers=sc_thin_layers,
                 small_niches=sc_small_niches, holes=sc_holes, variable_density=sc_variable_density,
                 gradient=sc_gradient)


def simulate(name, seed, noise=1.0):
    rng = np.random.default_rng(seed)
    xy, lab, cont = SCENARIOS[name](rng)
    K = len(np.unique(lab))
    prof = type_profiles(rng)
    if cont is not None:
        a, b = compositions(rng, 2, conc=0.5)
        comp = (1 - cont)[:, None] * a + cont[:, None] * b
    else:
        comp = compositions(rng, K)[lab]
    u = rng.uniform(size=len(xy))[:, None]
    ctype = (comp.cumsum(1) < u).sum(1).clip(0, N_TYPES - 1)
    X = prof[ctype] + rng.normal(0, noise, (len(xy), LAT))
    a = ad.AnnData(np.zeros((len(xy), 1), dtype=np.float32))
    a.obsm['spatial'] = xy; a.obsm['X_lat'] = X.astype(np.float32)
    a.obs['truth'] = lab; a.obs['ctype'] = ctype
    return a, cont


def boundary_mask(a, hops=2):
    # cells within `hops` graph steps of a cell with a different truth label
    A = a.obsp['spatial_connectivities'].tocsr(); A.data[:] = 1
    lab = a.obs['truth'].values
    diff = np.zeros(len(lab), bool)
    r, c = A.nonzero(); diff[r[lab[r] != lab[c]]] = True
    m = diff.copy()
    for _ in range(hops - 1):
        m = m | (A @ m.astype(float) > 0)
    return m


def matched_f1(truth, pred, cls):
    # Hungarian-match predicted clusters to truth, return F1 for a truth class
    ct = pd.crosstab(truth, pred)
    r, c = linear_sum_assignment(-ct.values)
    mp = {ct.columns[j]: ct.index[i] for i, j in zip(r, c)}
    mapped = np.array([mp.get(p, -1) for p in pred])
    return f1_score(truth == cls, mapped == cls)


def methods(a, n_layers):
    out = {}
    cc.gr.aggregate_neighbors(a, n_layers=n_layers, use_rep='X_lat', out_key='cc', sample_key=None)
    out['CellCharter'] = a.obsm['cc']
    for k in ['exponential', 'gaussian']:
        aggregate_neighbors_weighted(a, n_layers=n_layers, use_rep='X_lat', out_key=k,
                                     hop_decay=0.2, distance_kernel=k, normalize_weights=True)
        out[f'MANA-{k[:5]}'] = a.obsm[k]
    # hop decay applied the way it was presumably intended: scale each hop block AFTER normalisation
    blk = []
    for h in range(n_layers + 1):
        blk.append(a.obsm['gaussian'][:, h * LAT:(h + 1) * LAT] * (1.0 if h == 0 else 0.2 ** h / 0.2))
    out['MANA-gauss+decay'] = np.hstack(blk)
    return out


def run_one(name, seed, n_layers=3):
    a, cont = simulate(name, seed)
    sq.gr.spatial_neighbors(a, coord_type='generic', delaunay=True)
    cc.gr.remove_long_links(a)
    truth = a.obs['truth'].values; K = len(np.unique(truth))
    bmask = boundary_mask(a)
    rows = []
    for m, F in methods(a, n_layers).items():
        for clus in ['GMM', 'KMeans']:
            if clus == 'GMM':
                pred = GaussianMixture(K, covariance_type='full', n_init=3, random_state=seed).fit_predict(F)
            else:
                pred = KMeans(K, n_init=5, random_state=seed).fit_predict(F)
            row = dict(scenario=name, seed=seed, n_layers=n_layers, method=m, clusterer=clus,
                       ARI=adjusted_rand_score(truth, pred),
                       ARI_boundary=adjusted_rand_score(truth[bmask], pred[bmask]),
                       ARI_interior=adjusted_rand_score(truth[~bmask], pred[~bmask]))
            if name == 'small_niches':
                row['niche_F1'] = matched_f1(truth, pred, 2)
            if name == 'thin_layers':
                row['thin_F1'] = matched_f1(truth, pred, 1)
            rows.append(row)
        if cont is not None:  # how well do features linearly encode the gradient?
            r2 = cross_val_predict(RidgeCV(alphas=np.logspace(-2, 3, 8)), F, cont, cv=5)
            r2 = 1 - ((r2 - cont) ** 2).sum() / ((cont - cont.mean()) ** 2).sum()
            for r in rows[-2:]:
                r['gradient_R2'] = r2
    return rows


if __name__ == '__main__':
    import io, contextlib
    seeds = range(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
    layers = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [3]
    allrows = []
    for name in SCENARIOS:
        for L in layers:
            for s in seeds:
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    allrows += run_one(name, s, L)
            print('done', name, L, flush=True)
    pd.DataFrame(allrows).to_csv(sys.argv[3] if len(sys.argv) > 3 else 'results.csv', index=False)

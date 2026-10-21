"""Synthetic tissues with known domains and realistic within-domain structure.

The original simulator (synthetic_benchmark.py) draws cells independently within a domain, so
neighbour blocks shrink far more than in real tissue (self/3rd-ring variance ratio 12x vs 3.4x in
EAE) and effects such as feature scaling look much larger than they are. Here cell states vary
smoothly within domains (a random-Fourier-feature Gaussian field) and cell types cluster in
patches. Defaults (field_amp=1, field_len=6, patch_len=3, patch_strength=0.9) reproduce the EAE
ratios (1.96 / 2.61 / 3.45 vs EAE 2.23 / 2.92 / 3.44).

    python realistic_synthetic.py 10 results_realistic.csv              # all methods in methods.py
    python realistic_synthetic.py 10 out.csv "CellCharter,MANA shells"   # a subset
"""
import os, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, anndata as ad
from multiprocessing import Pool
from scipy.spatial import cKDTree
from scipy.stats import norm
from sklearn.metrics import adjusted_rand_score
import synthetic_benchmark as bench
from methods import METHODS, build_graph, cluster


def sc_lesion_rim(rng):
    xy = bench.jitter_grid(rng, 1.2)
    cen = rng.uniform(15, 85, (5, 2)); rad = rng.uniform(5, 10, 5)
    dd = (np.linalg.norm(xy[:, None] - cen[None], axis=2) - rad[None]).min(1)
    lab = np.digitize(dd, [0, 2.5, 10])  # core, rim (~2 cells), periplaque, normal
    return xy, lab, lab, np.eye(4)


def sc_patchy(rng):
    # left domain: mosaic of two compositions; middle: their 50/50 mix; right: another composition
    xy = bench.jitter_grid(rng, 1.2)
    lab = np.digitize(xy[:, 0], [40, 70])
    pseeds = rng.uniform(0, 40, (45, 2)) * [1, 2.5]
    patch = cKDTree(pseeds).query(xy)[1] % 2
    sub = np.where(lab == 0, patch, lab + 1)
    mix = np.array([[1, 0, 0], [0, 1, 0], [.5, .5, 0], [0, 0, 1]])
    return xy, lab, sub, mix


NEW = dict(lesion_rim=sc_lesion_rim, patchy=sc_patchy)
SCENARIOS = ['lesion_rim', 'patchy', 'holes', 'thin_layers', 'uniform', 'density_contrast',
             'small_niches', 'variable_density']


def smooth_field(rng, xy, length, n_out, M=300):
    """Gaussian-process-like field with unit variance per output (random Fourier features)."""
    W = rng.normal(0, 1 / length, (2, M)); b = rng.uniform(0, 2 * np.pi, M)
    return (np.sqrt(2 / M) * np.cos(xy @ W + b)) @ rng.normal(0, 1, (M, n_out))


def simulate(name, seed, noise=1.0, field_amp=1.0, field_len=6.0, patch_len=3.0, patch_strength=0.9):
    rng = np.random.default_rng(seed)
    if name in NEW:
        xy, lab, sub, mix = NEW[name](rng)
    else:
        xy, lab, _ = bench.SCENARIOS[name](rng); sub = lab; mix = np.eye(len(np.unique(lab)))
    prof = bench.type_profiles(rng)
    comp = (mix @ bench.compositions(rng, mix.shape[1], conc=0.5 if name in NEW else 0.6))[sub]
    rng2 = np.random.default_rng(seed + 10_000)
    z = patch_strength * smooth_field(rng2, xy, patch_len, 1)[:, 0] + np.sqrt(1 - patch_strength ** 2) * rng2.normal(size=len(xy))
    ctype = (comp.cumsum(1) < norm.cdf(z)[:, None]).sum(1).clip(0, bench.N_TYPES - 1)
    X = prof[ctype] + field_amp * smooth_field(rng2, xy, field_len, bench.LAT) + rng.normal(0, noise, (len(xy), bench.LAT))
    a = ad.AnnData(np.zeros((len(xy), 1), dtype=np.float32))
    a.obsm['spatial'] = xy; a.obsm['X_lat'] = X.astype(np.float32); a.obs['truth'] = lab
    return a


def run_one(args):
    name, seed, methods = args
    a = simulate(name, seed)
    build_graph(a)
    truth = a.obs['truth'].values; K = len(np.unique(truth))
    return [dict(scenario=name, seed=seed, method=m,
                 ARI=adjusted_rand_score(truth, cluster(METHODS[m](a, 'X_lat'), K, seed, threads=1)))
            for m in methods]


if __name__ == '__main__':
    n_seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    out = sys.argv[2] if len(sys.argv) > 2 else 'results_realistic.csv'
    methods = sys.argv[3].split(',') if len(sys.argv) > 3 else list(METHODS)
    with Pool(int(os.environ.get('N_JOBS', 8))) as p:
        rows = sum(p.map(run_one, [(s, i, methods) for s in SCENARIOS for i in range(n_seeds)], chunksize=1), [])
    d = pd.DataFrame(rows); d.to_csv(out, index=False)
    print(d.pivot_table(index='scenario', columns='method', values='ARI').round(3).to_string())
    print(d.groupby('method').ARI.mean().round(3).to_string())

"""Allen MERFISH whole mouse brain (MERFISH-C57BL6J-638850), scored against CCF atlas labels.

Labels come from registering each section to the Allen CCFv3. Fine labels (cortical layers,
CA1/CA3/DG, thalamic nuclei) follow the deformed atlas; major-region boundaries were partly fitted
to cell-type-derived region maps (Yao et al. 2023, Methods), so score mainly on the fine labels.

    python allen_merfish.py download DIR          # ~1.6 GB: 6 sections + coordinates/labels
    python allen_merfish.py prep DIR              # PCA, graph per section, features per method
    python allen_merfish.py bench DIR [methods]   # cc.tl.Cluster: cortex-only K=6 (10 seeds),
                                                  # whole sections K=20 (5 seeds)
"""
import os, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd

BUCKET = 'https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com'
SECTIONS = [30, 36, 38, 40, 44, 49]  # coronal, posterior hippocampus -> striatum
FILES = {
    'recon.csv': 'metadata/MERFISH-C57BL6J-638850-CCF/20231215/reconstructed_coordinates.csv',
    'parcellation.csv': 'metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_acronym.csv',
    **{f's{n}.h5ad': f'expression_matrices/MERFISH-C57BL6J-638850-sections/20240330/C57BL6J-638850.{n}-raw.h5ad'
       for n in SECTIONS},
}


def download(d):
    import urllib.request
    os.makedirs(d, exist_ok=True)
    for name, path in FILES.items():
        if not os.path.exists(os.path.join(d, name)):
            print('downloading', name, flush=True)
            urllib.request.urlretrieve(f'{BUCKET}/{path}', os.path.join(d, name))


def prep(d):
    import anndata as ad, scanpy as sc
    from methods import METHODS, build_graph
    parts = []
    for n in SECTIONS:
        b = ad.read_h5ad(os.path.join(d, f's{n}.h5ad'))
        b.obs_names = b.obs_names.astype(str); b.var_names = b.var['gene_symbol'].astype(str).values
        parts.append(b)
    a = ad.concat(parts); a.var_names_make_unique()
    rec = pd.read_csv(os.path.join(d, 'recon.csv'), usecols=['cell_label', 'x', 'y', 'parcellation_index'],
                      index_col='cell_label', dtype={'cell_label': str})
    rec = rec[rec.index.isin(a.obs_names)]
    a = a[rec.index].copy()
    lab = pd.read_csv(os.path.join(d, 'parcellation.csv'), index_col='parcellation_index').loc[rec['parcellation_index'].values]
    obs = pd.DataFrame(index=a.obs_names)
    obs['sample_id'] = a.obs['brain_section_label'].astype(str).values
    obs['division'] = lab['division'].values; obs['structure'] = lab['structure'].values
    obs['layer'] = np.where(lab['division'].values == 'Isocortex',
                            lab['substructure'].str.extract(r'(1|2/3|4|5|6a|6b)$')[0].values, 'None')
    obs['region'] = np.where(lab['category'].values == 'grey', lab['division'].values, lab['category'].values)
    sc.pp.normalize_total(a); sc.pp.log1p(a); sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=30)
    X = a.obsm['X_pca'].astype(np.float32); xy = rec[['x', 'y']].values * 1000  # mm -> um
    obs.astype(str).to_csv(os.path.join(d, 'obs.csv'))
    feats = {m: None for m in METHODS}
    for s, idx in obs.groupby('sample_id').indices.items():
        b = ad.AnnData(np.zeros((len(idx), 1), np.float32), obsm={'spatial': xy[idx], 'X_pca': X[idx]})
        build_graph(b)
        for m, f in METHODS.items():
            F = f(b, 'X_pca')
            if feats[m] is None:
                feats[m] = np.zeros((len(obs), F.shape[1]), np.float32)
            feats[m][idx] = F
        print('features done', s, flush=True)
    for m, F in feats.items():
        np.save(os.path.join(d, f'feat_{m.replace(" ", "_")}.npy'), F)


def _score(job):
    from sklearn.metrics import adjusted_rand_score as ARI, normalized_mutual_info_score as NMI
    from methods import cluster
    d, kind, m, K, seed = job
    obs = pd.read_csv(os.path.join(d, 'obs.csv'), index_col=0, keep_default_na=False)
    X = np.load(os.path.join(d, f'feat_{m.replace(" ", "_")}.npy'), mmap_mode='r')
    ctx = (obs.layer != 'None').values
    r = dict(kind=kind, method=m, K=K, seed=seed)
    if kind == 'cortex':
        lab = cluster(np.asarray(X[ctx]), K, seed, threads=3)
        r.update(ARI_layer=ARI(obs.layer.values[ctx], lab), NMI_layer=NMI(obs.layer.values[ctx], lab))
    else:
        lab = cluster(np.asarray(X), K, seed, threads=3)
        for nm, mask, col in [('layer', ctx, 'layer'),
                              ('hpf', obs.structure.isin(['CA1', 'CA2', 'CA3', 'DG']).values, 'structure'),
                              ('th', (obs.division == 'TH').values, 'structure'),
                              ('region', (obs.region != 'unassigned').values, 'region')]:
            r[f'ARI_{nm}'] = ARI(obs[col].values[mask], lab[mask]); r[f'NMI_{nm}'] = NMI(obs[col].values[mask], lab[mask])
    print(r, flush=True)
    return r


def bench(d, methods):
    from multiprocessing import Pool
    jobs = [(d, 'cortex', m, 6, s) for m in methods for s in range(10)] + \
           [(d, 'whole', m, 20, s) for m in methods for s in range(5)]
    with Pool(int(os.environ.get('N_JOBS', 5))) as p:
        res = pd.DataFrame(p.map(_score, jobs, chunksize=1))
    res.to_csv(os.path.join(d, 'allen_results.csv'), index=False)
    print(res.groupby(['kind', 'K', 'method']).mean(numeric_only=True).drop(columns='seed').round(3).to_string())


if __name__ == '__main__':
    step, d = sys.argv[1], sys.argv[2]
    if step == 'download':
        download(d)
    elif step == 'prep':
        prep(d)
    else:
        from methods import METHODS
        bench(d, sys.argv[3].split(',') if len(sys.argv) > 3 else list(METHODS))

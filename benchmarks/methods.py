"""Feature builders compared in the benchmarks. Add a new idea here and every benchmark picks it up.

Each builder takes one sample's AnnData (spatial graph already built with squidpy Delaunay +
cc.gr.remove_long_links) and the obsm key of the cell embedding, and returns a float32 array.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import io, contextlib
import numpy as np
from utils.aggregate_neighbors_weighted import aggregate_neighbors_weighted, aggregate_neighbors_shells


def _quiet(f, *args, **kw):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        f(*args, **kw)


def hops(a, rep, n_var=0, n_layers=3):
    # distance_kernel='none' == cc.gr.aggregate_neighbors (and avoids its scipy>=1.17 crash on >46K cells)
    _quiet(aggregate_neighbors_weighted, a, n_layers=n_layers, use_rep=rep, out_key='_f',
           distance_kernel='none', n_var_layers=n_var)
    return np.asarray(a.obsm.pop('_f'), dtype=np.float32)


def shells(a, rep, n_var=0, n_shells=3):
    _quiet(aggregate_neighbors_shells, a, n_shells=n_shells, use_rep=rep, out_key='_f', n_var_shells=n_var)
    return np.asarray(a.obsm.pop('_f'), dtype=np.float32)


METHODS = {
    'CellCharter': lambda a, rep: hops(a, rep),
    'CellCharter mean+var': lambda a, rep: hops(a, rep, n_var=3),   # == aggregations=['mean','var']
    'MANA shells': lambda a, rep: shells(a, rep),
    'MANA shells+var': lambda a, rep: shells(a, rep, n_var=3),
}


def build_graph(a):
    import squidpy as sq, cellcharter as cc
    _quiet(sq.gr.spatial_neighbors, a, coord_type='generic', delaunay=True)
    _quiet(cc.gr.remove_long_links, a)


def cluster(X, K, seed, threads=2):
    """cc.tl.Cluster (full-covariance GMM). Not reproducible across runs: always use several seeds."""
    import anndata as ad, cellcharter as cc, torch
    torch.set_num_threads(threads)
    b = ad.AnnData(np.zeros((len(X), 1), np.float32)); b.obsm['f'] = np.asarray(X, dtype=np.float32)
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        m = cc.tl.Cluster(n_clusters=K, random_state=seed,
                          trainer_params={'accelerator': 'cpu', 'enable_progress_bar': False})
        m.fit(b, use_rep='f')
        return np.asarray(m.predict(b, use_rep='f')).astype(int)


def smooth_labels(lab, A, iters=2):
    """Majority vote over each cell and its graph neighbours (A: spatial_connectivities), `iters` rounds.
    Annotation-free post-processing; ties go to the lowest label."""
    from scipy.sparse import csr_matrix, eye
    lab = np.asarray(lab); n = len(lab); K = lab.max() + 1
    Ahat = (A.astype(bool).astype(np.float32) + eye(n, format='csr')).tocsr()
    for _ in range(iters):
        onehot = csr_matrix((np.ones(n), (np.arange(n), lab)), shape=(n, K))
        lab = np.asarray((Ahat @ onehot).argmax(1)).ravel()
    return lab


def consensus_labels(runs):
    """Align each run's labels to runs[0] (Hungarian on the contingency table), then majority vote."""
    import pandas as pd
    from scipy.optimize import linear_sum_assignment
    ref = runs[0]; K = max(r.max() for r in runs) + 1
    votes = np.zeros((len(ref), K), np.float32)
    for r in runs:
        ct = pd.crosstab(r, ref)
        ri, ci = linear_sum_assignment(-ct.values)
        mp = np.arange(K); mp[ct.index[ri]] = ct.columns[ci]
        votes[np.arange(len(r)), mp[r]] += 1
    return votes.argmax(1)

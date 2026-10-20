"""Candidate improvements to MANA's within-neighbourhood weighting.
All return [self | block_1 | ... | block_n] like CellCharter, so dimensions match (n_blocks = n_layers)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np, pandas as pd
from scipy.sparse import csr_matrix, diags
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree
from utils.aggregate_neighbors_weighted import compute_hop_matrices


def _rownorm(W):
    s = np.asarray(W.sum(1)).ravel(); s[s == 0] = 1
    return diags(1 / s) @ W


def _binary_graph(a):
    A = a.obsp['spatial_connectivities'].tocsr().copy(); A.data[:] = 1
    return A


def _gauss(d, s):
    return np.exp(-d ** 2 / (2 * s ** 2))


def _row_median(rows, vals, n):
    med = pd.Series(vals).groupby(rows).transform('median').values
    return med


def geodesic_within(a, n_layers):
    """Shortest-path length along graph edges, for pairs within n_layers hops. Returns csr."""
    D = a.obsp['spatial_distances'].tocsr()
    A = _binary_graph(a)
    reach = A.copy()
    F = A.copy()
    for _ in range(n_layers - 1):
        F = F @ A; F.data[:] = 1; reach = reach + F; reach.data[:] = 1
    reach.setdiag(0); reach.eliminate_zeros()
    lim = np.percentile(D.data, 99) * n_layers
    n = D.shape[0]; out_r, out_c, out_d = [], [], []
    for start in range(0, n, 1500):
        idx = np.arange(start, min(n, start + 1500))
        G = dijkstra(D, indices=idx, limit=lim)
        sub = reach[idx]
        r, c = sub.nonzero()
        g = G[r, c]
        ok = np.isfinite(g)
        out_r.append(idx[r[ok]]); out_c.append(c[ok]); out_d.append(g[ok])
    return csr_matrix((np.concatenate(out_d), (np.concatenate(out_r), np.concatenate(out_c))), shape=(n, n))


def features(a, n_layers, kind, use_rep='X_lat'):
    X = np.asarray(a.obsm[use_rep]); xy = np.asarray(a.obsm['spatial'])
    A = _binary_graph(a)
    blocks = [X]
    if kind in ('blend', 'adaptive', 'geo_hop'):
        hops = compute_hop_matrices(A, n_layers)
        G = geodesic_within(a, n_layers) if kind == 'geo_hop' else None
        for H in hops:
            r, c = H.nonzero()
            if kind == 'geo_hop':
                d = np.asarray(G[r, c]).ravel()
                d[d == 0] = np.sqrt(((xy[r] - xy[c]) ** 2).sum(1))[d == 0]
            else:
                d = np.sqrt(((xy[r] - xy[c]) ** 2).sum(1))
            if kind == 'adaptive':
                s = _row_median(r, d, len(X))
            else:
                s = np.median(d) if len(d) else 1.0
            W = _rownorm(csr_matrix((_gauss(d, s), (r, c)), shape=H.shape))
            if kind == 'blend':
                W = 0.5 * W + 0.5 * _rownorm(H.astype(float))
            blocks.append(W @ X)
    elif kind in ('geo_multiscale', 'phys_multiscale', 'geo_shells', 'phys_shells'):
        e = np.median(a.obsp['spatial_distances'].tocsr().data)
        if kind.startswith('geo'):
            G = geodesic_within(a, n_layers)
            r, c = G.nonzero(); d = np.asarray(G[r, c]).ravel()
        else:
            P = cKDTree(xy).sparse_distance_matrix(cKDTree(xy), 3 * n_layers * e, output_type='coo_matrix')
            keep = P.row != P.col
            r, c, d = P.row[keep], P.col[keep], P.data[keep]
        for k in range(1, n_layers + 1):  # bandwidths e, 2e, 3e  ~ hop 1, 2, 3
            # shells: soft ring centred at k*e (continuous analogue of hop k); multiscale: nested blob
            w = _gauss(d - k * e, 0.5 * e) if kind.endswith('shells') else _gauss(d, k * e)
            W = _rownorm(csr_matrix((w, (r, c)), shape=(len(X), len(X))))
            blocks.append(W @ X)
    elif kind == 'phys_shells_adapt':
        # shells in units of each cell's own local spacing e_i (median Delaunay edge, smoothed over 1-hop)
        D = a.obsp['spatial_distances'].tocsr()
        ei = np.array([np.median(D.data[D.indptr[i]:D.indptr[i + 1]]) if D.indptr[i + 1] > D.indptr[i] else np.nan
                       for i in range(D.shape[0])])
        ei = np.where(np.isnan(ei), np.nanmedian(ei), ei)
        ei = 0.5 * ei + 0.5 * (_rownorm(A) @ ei)
        tree = cKDTree(xy)
        nb = tree.query_ball_point(xy, (n_layers + 1) * ei)
        r = np.repeat(np.arange(len(X)), [len(x) for x in nb]); c = np.concatenate(nb)
        keep = r != c; r, c = r[keep], c[keep]
        d = np.sqrt(((xy[r] - xy[c]) ** 2).sum(1)) / ei[r]  # distance in local-spacing units
        for k in range(1, n_layers + 1):
            W = _rownorm(csr_matrix((_gauss(d - k, 0.5), (r, c)), shape=(len(X), len(X))))
            blocks.append(W @ X)
    elif kind.startswith('alpha'):
        # partial density adaptation: e_i = e_global^(1-alpha) * e_local^alpha, e_local from 15-NN distance
        alpha = float(kind[5:])
        tree = cKDTree(xy)
        knn_d, _ = tree.query(xy, k=16)
        loc = knn_d[:, 1:].mean(1)
        e = np.median(a.obsp['spatial_distances'].tocsr().data)
        loc = loc / np.median(loc) * e  # same median as the global spacing
        loc = _rownorm(A) @ loc * 0.5 + loc * 0.5
        ei = e ** (1 - alpha) * loc ** alpha
        nb = tree.query_ball_point(xy, (n_layers + 1.5) * ei)
        r = np.repeat(np.arange(len(X)), [len(x) for x in nb]); c = np.concatenate(nb)
        keep = r != c; r, c = r[keep], c[keep]
        d = np.sqrt(((xy[r] - xy[c]) ** 2).sum(1)) / ei[r]
        for k in range(1, n_layers + 1):
            W = _rownorm(csr_matrix((_gauss(d - k, 0.5), (r, c)), shape=(len(X), len(X))))
            blocks.append(W @ X)
    elif kind.startswith('repo_shells'):
        from utils.aggregate_neighbors_weighted import aggregate_neighbors_shells
        aggregate_neighbors_shells(a, n_shells=n_layers, use_rep=use_rep, out_key='_sh_tmp',
                                   adaptive=float(kind.split('_')[-1]))
        return a.obsm.pop('_sh_tmp')
    elif kind == 'hybrid':
        # CellCharter hop rings (density-adaptive by construction) + physical shells (geometry-aware)
        import cellcharter as cc
        cc.gr.aggregate_neighbors(a, n_layers=n_layers, use_rep=use_rep, out_key='_cc_tmp')
        hop = a.obsm.pop('_cc_tmp')[:, X.shape[1]:]
        return np.hstack([features(a, n_layers, 'phys_shells', use_rep), hop])
    else:
        raise ValueError(kind)
    return np.hstack(blocks)

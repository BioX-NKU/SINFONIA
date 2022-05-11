#!/usr/bin/env python

import numpy as np
import scipy.sparse as sp
from numba import jit
from sklearn.neighbors import kneighbors_graph
from .utils import *

def spatially_variable_genes(adata, mode='moran_geary', n_top_genes=2000, subset=False, inplace=True, 
    coordinate_key='spatial', n_neighbors=30, layer=None, seed=None):
    """
    Identify spatially variable genes by SINFONIA.
    
    Mode
    -----
        * 1. `moran_geary` identify spatially variable genes based on Moran’s I and rescaled Geary’s C scores.
        * 2. `moran` identify spatially variable genes based on Moran’s I scores.
        * 3. `geary` identify spatially variable genes based on rescaled Geary’s C scores.

    Parameters
    ----------
    adata
        AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    mode
        Mode for identifying spatially variable genes. By default, mode='moran_geary'.
    n_top_genes
        Number of spatially variable genes to keep. By default, n_top_genes=2000.
    subset
        Inplace subset to spatially variable genes if `True` otherwise merely indicate spatially variable genes. By default, subset=False.
    inplace
        Whether to place calculated metrics in `.var`. By default, inplace=True.
    coordinate_key
        Key in `adata.obsm` where spatial coordinates are stored. By default, coordinate_key='spatial'.
    n_neighbors
        Number of nearest neighbors used to construct the spatial neighbor graph. By default, n_neighbors=30.
    layer
        If provided, use `adata.layers[layer]` for expression values instead of `adata.X`. Note that SINFONIA expects logarithmized data. By default, layer=None.
    seed
        Random seed for reproducibility. By default, seed=None.
        
    Returns
    -------
    adata
        AnnData object with identified spatially variable genes and additional fields:
        - `adata.var['spatially_variable']` - boolean indicator of spatially variable genes.
        - `adata.uns['svg']` - mode used to identify spatially variable genes if `inplace=True`.
        - `adata.var['moranI']` - Moran’s I scores if `inplace=True` and `mode='moran_geary'` or `mode='moran'`.
        - `adata.var['gearyC']` - rescaled Geary’s C scores if `inplace=True` and `mode='moran_geary'` or `mode='geary'`.

    """

    if seed is not None: setup_seed(seed)
    sc.pp.filter_genes(adata, min_cells=1)

    @jit(nopython=True, parallel=True)
    def _moran(exp_mtx, data, row, col):
        N = exp_mtx.shape[0]
        sum_W = np.sum(data)
        moran_index = []
        for i in range(exp_mtx.shape[1]):
            exp_vec = exp_mtx[:, i]
            vec_mean = exp_vec.mean()
            C = N*(data*(exp_vec[row]-vec_mean)*(exp_vec[col]-vec_mean)).sum() / (sum_W*np.sum((exp_vec-vec_mean)**2))
            moran_index.append(C)
        return moran_index
    
    @jit(nopython=True, parallel=True)
    def _geary(exp_mtx, data, row, col):
        N = exp_mtx.shape[0]
        sum_W = np.sum(data)
        geary_index = []
        for i in range(exp_mtx.shape[1]):
            exp_vec = exp_mtx[:, i]
            C = (N-1)*(data*((exp_vec[row]-exp_vec[col])**2)).sum() / (2*sum_W*np.sum((exp_vec-exp_vec.mean())**2))
            geary_index.append(1-C)
        return geary_index
    
    @jit(nopython=True, parallel=True)
    def _moran_geary(exp_mtx, data, row, col):
        N = exp_mtx.shape[0]
        sum_W = np.sum(data)
        moran_index = []
        geary_index = []
        for i in range(exp_mtx.shape[1]):
            exp_vec = exp_mtx[:, i]
            vec_mean = exp_vec.mean()
            down_item = np.sum((exp_vec-vec_mean)**2)
            exp_row, exp_col = exp_vec[row], exp_vec[col]
            I = N*(data*(exp_row-vec_mean)*(exp_col-vec_mean)).sum() / (sum_W*down_item)
            C = (N-1)*(data*((exp_row-exp_col)**2)).sum() / (2*sum_W*down_item)
            moran_index.append(I)
            geary_index.append(1-C)
        return moran_index, geary_index
    
    W = kneighbors_graph(adata.obsm[coordinate_key], n_neighbors=n_neighbors,
                        include_self=False, metric="euclidean", mode="distance").tocoo()
    W.data = 1 - W.data/(W.max(axis=1).data.reshape(-1, 1).
                         repeat(n_neighbors, axis=1).reshape(-1) + 1e-5)
    X = adata.layers[layer] if layer is not None else adata.X
    exp_mtx = X.toarray() if sp.issparse(X) else X.copy()
    
    if mode=='moran_geary':
        moran_index, geary_index = _moran_geary(exp_mtx, W.data, W.row, W.col)
        genes_selected1 = np.array([False] * adata.shape[1])
        genes_selected1[np.argpartition(moran_index, len(moran_index)-n_top_genes)[-n_top_genes:]] = True
        genes_selected2 = np.array([False] * adata.shape[1])
        genes_selected2[np.argpartition(geary_index, len(geary_index)-n_top_genes)[-n_top_genes:]] = True    
        genes_selected = np.logical_or(genes_selected1, genes_selected2)
    elif mode=='moran':
        moran_index = _moran(exp_mtx, W.data, W.row, W.col)
        genes_selected = np.array([False] * adata.shape[1])
        genes_selected[np.argpartition(moran_index, len(moran_index)-n_top_genes)[-n_top_genes:]] = True
    elif mode=='geary':    
        geary_index = _geary(exp_mtx, W.data, W.row, W.col)
        genes_selected = np.array([False] * adata.shape[1])
        genes_selected[np.argpartition(geary_index, len(geary_index)-n_top_genes)[-n_top_genes:]] = True
    adata.var['spatially_variable'] = genes_selected
    
    if inplace:
        adata.uns['svg'] = {'mode': mode}
        if mode=='moran_geary':
            adata.var['moranI'] = moran_index
            adata.var['gearyC'] = geary_index
        elif mode=='moran':
            adata.var['moranI'] = moran_index
        elif mode=='geary':    
            adata.var['gearyC'] = geary_index
    if subset:
        adata = adata[:, genes_selected]

    return adata


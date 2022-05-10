#!/usr/bin/env python

import random
import numpy as np
import scanpy as sc

def setup_seed(seed):
    """
    Set random seed.
    
    Parameters
    ----------
    seed
        Number to be set as random seed for reproducibility.
        
    """
    np.random.seed(seed)
    random.seed(seed)

def get_N_clusters(adata, n_cluster, cluster_method='louvain', range_min=0, range_max=3, max_steps=30, tolerance=0):
    """
    Tune the resolution parameter in clustering to make the number of clusters and the specified number as close as possible.
    
    Parameters
    ----------
    adata
        AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    n_cluster
        Specified number of clusters.
    cluster_method
        Method (`louvain` or `leiden`) used for clustering. By default, cluster_method='louvain'.
    range_min
        Minimum clustering resolution for the binary search.
    range_max
        Maximum clustering resolution for the binary search.
    max_steps
        Maximum number of steps for the binary search.
    tolerance
        Tolerance of the difference between the number of clusters and the specified number.

    Returns
    -------
    adata
        AnnData object with clustering assignments in `adata.obs`:
        - `adata.obs['louvain']` - Louvain clustering assignments if `cluster_method='louvain'`.
        - `adata.obs['leiden']` - Leiden clustering assignments if `cluster_method='leiden'`.

    """ 
    this_step = 0
    this_min = float(range_min)
    this_max = float(range_max)
    while this_step < max_steps:
        this_resolution = this_min + ((this_max-this_min)/2)
        if cluster_method=='leiden':
            sc.tl.leiden(adata, resolution=this_resolution)
        if cluster_method=='louvain':
            sc.tl.louvain(adata, resolution=this_resolution)
        this_clusters = adata.obs[cluster_method].nunique()

        if this_clusters > n_cluster+tolerance:
            this_max = this_resolution
        elif this_clusters < n_cluster-tolerance:
            this_min = this_resolution
        else:
            print("Succeed to find %d clusters at resolution %.3f."%(n_cluster, this_resolution))
            return adata
        this_step += 1
    
    print('Cannot find the number of clusters.')
    return adata

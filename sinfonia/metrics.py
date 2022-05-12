#!/usr/bin/env python

import numpy as np
import pandas as pd
import sklearn
from sklearn.svm import SVC
from sklearn.model_selection import cross_validate
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import homogeneity_score
from sklearn.metrics.cluster import normalized_mutual_info_score


def clustering_metrics(adata, target, pred):   
    """
    Evaluate clustering performance.
    
    Parameters
    ----------
    adata
        AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
    target
        Key in `adata.obs` where ground-truth spatial domain labels are stored.
    pred
        Key in `adata.obs` where clustering assignments are stored.
        
    Returns
    -------
    ami
        Adjusted mutual information score.
    ari
        Adjusted Rand index score.
    homo
        Homogeneity score.
    nmi
        Normalized mutual information score.

    """ 
    ami  = adjusted_mutual_info_score(adata.obs[target], adata.obs[pred])
    ari  = adjusted_rand_score(adata.obs[target], adata.obs[pred])
    homo = homogeneity_score(adata.obs[target], adata.obs[pred])
    nmi  = normalized_mutual_info_score(adata.obs[target], adata.obs[pred])

    return ami, ari, homo, nmi

def mean_average_precision(embedding: np.ndarray, labels: np.ndarray, n_neighbors: int=30, **kwargs) -> float:
    """
    Evaluate mean average precision.

    Parameters
    ----------
    embedding
        Low-dimensional representation of spots.
    labels
        Ground-truth domain labels.
    n_neighbors
        Number of nearest neighbors used to evaluate mean average precision. By default, n_neighbors=30.
    **kwargs
        Additional keyword arguments are passed to
        :class:`sklearn.neighbors.NearestNeighbors`
        
    Returns
    -------
    map
        Mean average precision.
        
    """
    def _average_precision(match: np.ndarray) -> float:
        if np.any(match):
            cummean = np.cumsum(match) / (np.arange(match.size) + 1)
            return cummean[match].mean().item()
        return 0.0
    if isinstance(labels, pd.Series):
        labels = labels.astype(str).values
    knn = sklearn.neighbors.NearestNeighbors(
        n_neighbors=min(labels.shape[0], n_neighbors + 1), **kwargs
    ).fit(embedding)
    nni = knn.kneighbors(embedding, return_distance=False)
    match = np.equal(labels[nni[:, 1:]], np.expand_dims(labels, 1))
    
    return np.apply_along_axis(_average_precision, 1, match).mean().item()

def mean_cross_validation_accuracy(embedding, labels, Kfold=5):
    """
    Evaluate mean cross-validation accuracy.

    Parameters
    ----------
    embedding
        Low-dimensional representation of spots.
    labels
        Ground-truth domain labels.
    Kfold
        Number of folds for cross-validation. By default, Kfold=5.

    Returns
    -------
    MCVA
        Mean cross-validation accuracy.

    """

    if isinstance(labels, pd.Series):
        labels = labels.astype(str).values
    target_unique = np.unique(labels).reshape(1, -1)
    target_onehot = (labels.reshape(-1, 1)==target_unique).astype(int)
    labels = target_onehot.argmax(-1)
    svc = SVC()
    cv_results = cross_validate(svc, embedding, labels, scoring="accuracy", cv=Kfold, n_jobs=Kfold)
    
    return np.mean(cv_results["test_score"])
    
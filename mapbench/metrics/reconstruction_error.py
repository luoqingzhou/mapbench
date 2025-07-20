import numpy as np
import pandas as pd
import scipy
from typing import Union
from anndata import AnnData
from scvi.model import SCVI
from sklearn.metrics.pairwise import cosine_distances

def reconstruction_error(
    adata: AnnData,
    model_path: str,
    ref_query_key: str = "ref_query",
    query_key: str = "query",
    n_samples: int = 50,
    seed: int = 42,
    use_cosine: bool = True,
    batch_key: str = "sample_id",
    save_to_adata: bool = True,
) -> Union[AnnData, np.ndarray]:
    """
    Compute reconstruction error on query cells using a pretrained SCVI model.
    """
    import scvi
    scvi.settings.seed = seed

    # Load model and bind query data
    model = SCVI.load(model_path, adata=adata)
    query_adata = adata[adata.obs[ref_query_key] == query_key].copy()
    model = model.load_query_data(query_adata)

    # Normalization if needed
    if "log1p" not in query_adata.uns:
        import scanpy as sc
        sc.pp.normalize_total(query_adata, target_sum=10000, layer='raw')
        sc.pp.log1p(query_adata)

    # Get original data
    X_true = query_adata[:, model.adata.var_names].layers['raw'].copy()
    posterior_samples = model.posterior_predictive_sample(
        indices=np.arange(query_adata.n_obs), n_samples=n_samples
    )

    normalized_samples = np.zeros_like(posterior_samples)
    for s in range(posterior_samples.shape[2]):
        normalized_samples[:, :, s] = (
            (posterior_samples[:, :, s].T / posterior_samples[:, :, s].sum(axis=1) * 10000).T
        )

    log_samples = np.log1p(normalized_samples)
    X_pred = log_samples.mean(axis=2)

    if scipy.sparse.issparse(X_true):
        X_true = X_true.toarray()
    X_pred = X_pred.astype(X_true.dtype)

    if use_cosine:
        distances = cosine_distances(X_true, X_pred)
        errors = np.diag(distances)
    else:
        errors = np.mean((X_true - X_pred) ** 2, axis=1)

    if save_to_adata:
        query_adata.obs['reconstruction_error'] = errors
        query_adata.uns['reconstruction_model'] = model_path
        return query_adata
    else:
        return errors
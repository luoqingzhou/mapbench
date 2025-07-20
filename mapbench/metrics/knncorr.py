import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors
import anndata as ad
import scipy

def knncorr(
    adata: ad.AnnData,
    ref_query_key: str = "ref_query",
    query_key: str = "query",
    embedding: str = "X",
    k: int = 100,
    distance: str = "euclidean"
) -> np.ndarray:
    """
    Calculates k-NN correlation within the query cells only.
    
    Parameters:
    - adata: AnnData object that contains both reference and query cells.
    - ref_query_key: column in `.obs` to distinguish reference/query cells.
    - query_key: value in `ref_query_key` corresponding to query cells.
    - embedding: name of the .obsm slot with reference mapping coordinates (e.g., "X").
    - k: number of nearest neighbors to use.
    - distance: either 'euclidean' or 'cosine'.
    
    Returns:
    - A 1D numpy array of Spearman correlations for each query cell.
    """

    # Get query cells
    if ref_query_key not in adata.obs.columns:
        raise ValueError(f"Missing ref/query indicator column `{ref_query_key}` in adata.obs")
    query_mask = adata.obs[ref_query_key] == query_key
    query_adata = adata[query_mask]

    if "X_pca" not in adata.obsm:
        raise ValueError("Missing `.obsm['X_pca']` — this PCA representation is required for calculating ground-truth embedding.")
    
    X_pca = query_adata.obsm["X_pca"].toarray() if scipy.sparse.issparse(query_adata.obsm["X_pca"]) else query_adata.obsm["X_pca"]
    if embedding == "X":
        X_map = query_adata.X.toarray() if scipy.sparse.issparse(query_adata.X) else query_adata.X
    else:
        if embedding not in adata.obsm:
            raise ValueError(f"Missing `.obsm['{embedding}']` — this reference mapping embedding is required.")
        else:
            X_map = query_adata.obsm[embedding].toarray() if scipy.sparse.issparse(query_adata.obsm[embedding]) else query_adata.obsm[embedding]


    if X_pca.shape != X_map.shape:
        raise ValueError(f"Shape mismatch between X_pca {X_pca.shape} and embedding {X_map.shape}")

    n_cells = X_pca.shape[0]
    if n_cells <= k:
        k = n_cells - 1
        print(f"Warning: too few query cells ({n_cells}). Using k={k}")

    # Normalize if using cosine
    if distance == "cosine":
        def normalize(x): return x / np.linalg.norm(x, axis=1, keepdims=True)
        X_pca = normalize(X_pca)
        X_map = normalize(X_map)
    elif distance != "euclidean":
        raise ValueError("Unsupported distance metric: choose from 'euclidean' or 'cosine'.")

    # Fit NN on PCA space
    nbrs = NearestNeighbors(n_neighbors=k + 1, metric=distance).fit(X_pca)
    knn_indices = nbrs.kneighbors(X_pca, return_distance=False)[:, 1:]  # exclude self

    corrs = np.zeros(n_cells)
    for i in range(n_cells):
        idx = knn_indices[i]

        anchor_pca = np.tile(X_pca[i], (k, 1))
        anchor_map = np.tile(X_map[i], (k, 1))
        neighbors_pca = X_pca[idx]
        neighbors_map = X_map[idx]
        
        dist_pca = np.linalg.norm(anchor_pca - neighbors_pca, axis=1)
        dist_map = np.linalg.norm(anchor_map - neighbors_map, axis=1)

        rho, _ = spearmanr(dist_pca, dist_map)
        corrs[i] = rho if not np.isnan(rho) else 0.0

    return corrs

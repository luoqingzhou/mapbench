import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsTransformer

def uncertainty(
    adata,
    embedding="X",
    label_key="cell_type",
    ref_query_key="ref_query",
    ref_key="ref",
    query_key="query",
    n_neighbors=50,
    threshold=1.0,
    pred_unknown=False,
):
    """
    Calculate uncertainty metric using weighted KNN between ref and query cells.

    Parameters
    ----------
    adata : AnnData
        Combined AnnData object containing both reference and query cells.
    embedding : str
        Name of obsm layer or 'X' to use for embeddings.
    label_key : str
        Column in .obs representing cell labels.
    ref_query_key : str
        Column in .obs indicating reference/query category.
    ref_key : str
        Value in ref_query_key identifying reference cells.
    query_key : str
        Value in ref_query_key identifying query cells.
    n_neighbors : int
        Number of neighbors to consider.
    threshold : float
        Uncertainty threshold for marking as unknown.
    pred_unknown : bool
        Whether to assign 'Unknown' if below threshold.

    Returns
    -------
    pred_labels : pd.Series
        Predicted labels for query cells.
    uncertainties : pd.Series
        Uncertainty scores for query cells.
    """
    ref_adata = adata[adata.obs[ref_query_key] == ref_key]
    query_adata = adata[adata.obs[ref_query_key] == query_key]

    if embedding == "X":
        ref_emb = ref_adata.X
        query_emb = query_adata.X
    else:
        ref_emb = ref_adata.obsm[embedding]
        query_emb = query_adata.obsm[embedding]

    y_train_labels = ref_adata.obs[label_key].values

    # Train weighted KNN on reference
    knn = KNeighborsTransformer(
        n_neighbors=n_neighbors,
        mode="distance",
        algorithm="brute",
        metric="euclidean",
        n_jobs=-1,
    )
    knn.fit(ref_emb)

    top_k_distances, top_k_indices = knn.kneighbors(query_emb)

    stds = np.std(top_k_distances, axis=1)
    stds = (2.0 / stds) ** 2
    stds = stds.reshape(-1, 1)

    weights = np.exp(-top_k_distances / stds)
    weights /= np.sum(weights, axis=1, keepdims=True)

    pred_labels = []
    uncertainties = []

    for i in range(len(query_adata)):
        neighbor_labels = y_train_labels[top_k_indices[i]]
        unique_labels = np.unique(neighbor_labels)

        best_label = None
        best_prob = 0.0

        for candidate_label in unique_labels:
            candidate_prob = weights[i, neighbor_labels == candidate_label].sum()
            if candidate_prob > best_prob:
                best_prob = candidate_prob
                best_label = candidate_label

        if pred_unknown and best_prob < threshold:
            pred_labels.append("Unknown")
        else:
            pred_labels.append(best_label)

        uncertainties.append(max(1 - best_prob, 0))

    pred_labels = pd.Series(pred_labels, index=query_adata.obs_names, name="pred_label")
    uncertainties = pd.Series(uncertainties, index=query_adata.obs_names, name="uncertainty")

    print("Uncertainty calculation finished!")

    # return pred_labels, uncertainties
    return uncertainties
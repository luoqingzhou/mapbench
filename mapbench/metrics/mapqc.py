from mapqc import run_mapqc

def MapQC(adata, embedding="X",ref_query_key="ref_query", ref_key="ref", query_key="query", batch_key="sample_id", study_key="dataset_id", **kwargs):
    """
    Run MapQC on the given AnnData object.

    Parameters:
    - adata: AnnData object containing the data.
    - embedding: Key for the embedding to use (default: "X").
    - label_key: Key for the labels in adata.obs (default: "cell_type").
    - ref_query_key: Key to distinguish reference and query cells (default: "ref_query").
    - ref_key: Value in ref_query_key for reference cells (default: "ref").
    - query_key: Value in ref_query_key for query cells (default: "query").
    - n_neighbors: Number of neighbors to consider for uncertainty calculation (default: 5).
    - kwargs: Additional keyword arguments for run_mapqc.

    Returns:
    - adata with MapQC results added to uns.
    """
    return run_mapqc(
        adata=adata, 
        adata_emb_loc=embedding, 
        ref_q_key=ref_query_key,
        q_cat=query_key,
        r_cat=ref_key,
        sample_key=batch_key,
        study_key=study_key,
        **kwargs
        )
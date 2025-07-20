import logging
import warnings

import milopy
import pandas as pd
import scanpy as sc
from anndata import AnnData

def run_milo(
    adata: AnnData,
    ref_query_key: str = "ref_query",
    ref_key: str = "ref",
    query_key: str = "query",
    batch_key: str = "sample_id",
    celltype_key: str = "cell_type",
    design: str = "~is_query",
):
    milopy.core.make_nhoods(adata, prop=0.1)
    milopy.core.count_nhoods(adata, sample_col=batch_key)
    milopy.utils.annotate_nhoods(adata[adata.obs[ref_query_key] == ref_key], celltype_key)
    adata.obs["is_query"] = adata.obs[ref_query_key] == query_key
    milopy.core.DA_nhoods(adata, design=design)


def DALogFC(
    adata: AnnData,
    embedding: str = "X",
    ref_query_key: str = "ref_query",
    ref_key: str = "ref",
    query_key: str = "query",
    batch_key: str = "sample_id",
    celltype_key: str = "cell_type",
    milo_design: str = "~is_query",
    **kwargs,
):
    # Make KNN graph for Milo neigbourhoods
    n_controls = adata[adata.obs[ref_query_key] == ref_key].obs[batch_key].unique().shape[0]
    n_querys = adata[adata.obs[ref_query_key] == query_key].obs[batch_key].unique().shape[0]
    #  Set max to 200 or memory explodes for large datasets
    k = min([(n_controls + n_querys) * 5, 200])
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata, use_rep=embedding, n_neighbors=k)
    run_milo(adata, ref_query_key, ref_key, query_key, batch_key, celltype_key, milo_design)

    sample_adata = adata.uns["nhood_adata"].T.copy()
    sample_adata.var["OOR_score"] = sample_adata.var["logFC"].copy()
    sample_adata.var["OOR_signif"] = (
        ((sample_adata.var["SpatialFDR"] < 0.1) & (sample_adata.var["logFC"] > 0)).astype(int).copy()
    )
    sample_adata.varm["groups"] = adata.obsm["nhoods"].T
    adata.uns["sample_adata"] = sample_adata.copy()
    return adata
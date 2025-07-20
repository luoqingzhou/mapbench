import scanpy as sc
import pandas as pd
import numpy as np
import mapqc
from mapbench.metrics import DALogFC
from mapbench.metrics import uncertainty
from mapbench.metrics import MapQC
from mapbench.evaluation import compute_identification_metrics
import os
pd.DataFrame.iteritems = pd.DataFrame.items


def run_evaluation(adata_path, celltype_key, batch_key, mapqc_config, output_dir="."):
    """
    Run evaluation using mapQC, DAlogFC, and Uncertainty on the given AnnData.
    Saves a single row of results to evaluation_result.csv

    Parameters:
    - adata: AnnData object (must contain .uns['oor_celltype'])
    """
    adata = sc.read_h5ad(adata_path)
    adata_name = os.path.basename(adata_path).replace(".h5ad", "")

    assert "oor_celltype" in adata.uns, "Missing `oor_celltype` in adata.uns"
    assert "ref_query" in adata.obs, "Missing `ref_query` in adata.obs"
    target_celltype = adata.uns["oor_celltype"]

    # ======== 添加 case/control 标签 ========
    adata.obs['case_control'] = 'not_query'  # 先初始化所有样本为 'not_query'（字符串类型）
    query_mask = adata.obs['ref_query'] == 'query'  # 筛选 query 样本
    adata.obs.loc[query_mask, 'case_control'] = 'control'  # query 样本默认设为 'control'

    # 目标细胞类型样本设为 'case'（仅在 query 样本中）
    target_mask = query_mask & (adata.obs[celltype_key] == target_celltype)
    adata.obs.loc[target_mask, 'case_control'] = 'case'

    # ======== 1. Uncertainty ========
    print("→ Computing uncertainty...")
    adata = uncertainty(adata, label_key=celltype_key)
    adata_query = adata[adata.obs["ref_query"] == "query"].copy()
    y_true = adata_query.obs[celltype_key] == target_celltype
    y_pred = adata.uns['uncertainty']
    if y_true.sum() == 0:
        print(f"Warning: No samples found for target cell type '{target_celltype}' in uncertainty scores.")
        unc_auprc = 0
    else:
        unc_auprc = compute_identification_metrics(y_true=y_true.values, y_score=y_pred.values)

    # ======== 2. DAlogFC ========
    print("→ Computing DAlogFC...")
    adata = DALogFC(adata, celltype_key=celltype_key, batch_key=batch_key)
    idx = adata.uns["nhood_adata"].obs["index_cell"]
    adata_subset = adata[adata.obs_names.isin(idx)].copy()
    y_true = adata_subset.obs[celltype_key] == target_celltype
    y_pred = adata.uns["nhood_adata"].obs["logFC"]
    if y_true.sum() == 0:
        print(f"Warning: No samples found for target cell type '{target_celltype}' in DAlogFC scores.")
        da_auprc = 0
    else:
        da_auprc = compute_identification_metrics(y_true=y_true.values, y_score=y_pred.values)

    n_nhoods, k_min, k_max, study_key = mapqc_config.get("n_nhoods", 300), mapqc_config.get("k_min", 1500), mapqc_config.get("k_max", 10000), mapqc_config.get("study_key", "dataset")
    # ======== 3. mapQC ========
    print("→ Computing mapQC...")
    MapQC(
        adata,
        embedding="X",
        ref_query_key="ref_query",
        ref_key="ref",
        query_key="query",
        study_key=study_key,
        batch_key=batch_key,
        n_nhoods=n_nhoods,
        k_min=k_min,
        k_max=k_max,
        seed=10,
        overwrite=True,
        return_nhood_info_df=False,
        return_sample_dists_to_ref_df=False,
    )

    mapqc_stats = mapqc.evaluate(
        adata,
        case_control_key="case_control",
        case_cats=["case"],
        control_cats=["control"]
    )
    adata_mapqc = adata[adata.obs["mapqc_score"].notna()].copy()
    y_true = adata_mapqc.obs[celltype_key] == target_celltype
    y_pred = adata_mapqc.obs["mapqc_score"]
    if y_true.sum() == 0:
        print(f"Warning: No samples found for target cell type '{target_celltype}' in mapQC scores.")
        mapqc_auprc = 0
    else:
        mapqc_auprc = compute_identification_metrics(y_true=y_true.values, y_score=y_pred.values)


    # ======== 写入 CSV（追加模式） ========
    row = pd.DataFrame([{
        "adata_name": adata_name,
        "Uncertainty": unc_auprc,
        "DAlogFC": da_auprc,
        "mapQC": mapqc_auprc
    }])

    output_path = os.path.join(output_dir, "evaluation_result.csv")
    if not os.path.exists(output_path):
        row.to_csv(output_path, index=False, sep='\t')
    else:
        row.to_csv(output_path, index=False, mode="a", header=False, sep='\t')
    return adata

import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler, RobustScaler
from scipy.spatial.distance import cdist

def dist_uncertainty(
    adata,
    embedding="X",  # 假设 scPoli 得到的 latent 存在这里
    label_key="cell_type",
    ref_query_key="ref_query",
    ref_key="ref",
    query_key="query",
    scale=True,
    return_adata=True
):
    """
    计算基于潜空间原型距离的不确定性指标 (Inspired by scPoli).
    
    指标逻辑: 
    Uncertainty = min(Distance(cell_i, Prototype_k))
    即细胞距离其最接近的已知类型中心的距离。
    """
    # 1. 提取参考集和查询集数据
    ref_mask = adata.obs[ref_query_key] == ref_key
    query_mask = adata.obs[ref_query_key] == query_key
    
    ref_adata = adata[ref_mask]
    query_adata = adata[query_mask]

    # 获取嵌入向量 (Embedding)
    if embedding == "X":
        ref_emb = ref_adata.X
        query_emb = query_adata.X
    else:
        ref_emb = ref_adata.obsm[embedding]
        query_emb = query_adata.obsm[embedding]

    # 确保是 dense array
    if hasattr(ref_emb, "toarray"): ref_emb = ref_emb.toarray()
    if hasattr(query_emb, "toarray"): query_emb = query_emb.toarray()

    # 2. 计算参考集中每个类别的原型 (Prototypes)
    # 这里的原型即为各细胞类型在潜空间的质心
    labels = ref_adata.obs[label_key]
    unique_labels = labels.unique()
    prototypes = []
    
    for label in unique_labels:
        mask = (labels == label).values
        centroid = ref_emb[mask].mean(axis=0)
        prototypes.append(centroid)
    
    prototypes = np.stack(prototypes) # 形状: [类别数, 维度]

    # 3. 计算 Query 细胞到所有原型的距离
    # distances 形状: [Query细胞数, 类别数]
    distances = cdist(query_emb, prototypes, metric='euclidean')

    # 4. 提取最小距离作为不确定性 (公式 15)
    # 并记录预测标签 (公式 14)
    min_dist = np.min(distances, axis=1)
    best_proto_idx = np.argmin(distances, axis=1)
    pred_labels = unique_labels[best_proto_idx]

    # 5. 标准化处理 (模仿 scPoli 源码逻辑)
    if scale:
        # 使用 RobustScaler 减少离群值影响，再映射到 0-1
        scaler_res = RobustScaler().fit_transform(min_dist.reshape(-1, 1))
        uncertainties = MinMaxScaler(feature_range=(0, 1)).fit_transform(scaler_res).flatten()
    else:
        uncertainties = min_dist

    # 封装结果
    res_uncertainty = pd.Series(uncertainties, index=query_adata.obs_names, name="dist_uncertainty")
    res_preds = pd.Series(pred_labels, index=query_adata.obs_names, name="dist_pred_label")

    print("Dist_uncertainty calculation finished!")

    # 写入 adata
    # 注意：由于只计算了 query 部分，其余部分填入 NaN
    adata.obs['dist_uncertainty'] = np.nan
    adata.obs.loc[query_mask, 'dist_uncertainty'] = res_uncertainty
    
    adata.obs['dist_pred_label'] = "Reference"
    adata.obs.loc[query_mask, 'dist_pred_label'] = res_preds

    adata.uns['dist_uncertainty'] = res_uncertainty

    if return_adata:
        return adata
    return res_uncertainty
import numpy as np
from sklearn.metrics import (
    roc_auc_score, 
    average_precision_score, 
    precision_recall_curve, 
    confusion_matrix
)

def compute_identification_metrics(y_true, y_score, target_fdr=0.1):
    """
    针对单细胞异常检测优化的评估函数。
    
    参数:
    y_true: 数组, 真实标签 (0: 正常/参考, 1: 异常/OOR)
    y_score: 数组, 算法输出的异常得分 (越高表示越异常)
    target_fdr: 目标假发现率 (默认 0.1，即控制 10% 的误报)
    
    返回:
    包含各项指标的平铺字典，便于直接转为 DataFrame 列。
    """
    
    # 强制转换为 numpy 数组并处理空值/异常
    y_true = np.array(y_true).ravel()
    y_score = np.array(y_score).ravel()
    
    # 边界情况处理：如果全是正常或全是异常，无法计算 AUC
    if len(np.unique(y_true)) < 2:
        return {
            "AUPRC": 0.0, "AUROC": 0.0,
            "TPR_at_10FDR": 0.0, "FPR_at_10FDR": 0.0,
            "Actual_Precision": 0.0, "Threshold_at_10FDR": 0.0
        }

    # 1. 整体性能指标 (整体区分度)
    auroc = roc_auc_score(y_true, y_score)
    auprc = average_precision_score(y_true, y_score)

    # 2. 寻找满足目标 FDR 的阈值
    # Precision = TP / (TP + FP)
    # FDR = FP / (TP + FP) = 1 - Precision
    # 目标 10% FDR 意味着需要 Precision >= 0.9
    precisions, recalls, thresholds = precision_recall_curve(y_true, y_score)
    
    # 寻找 Precision >= 1 - target_fdr 的所有阈值索引
    # 注意：precisions 最后一项是 1，recalls 最后一项是 0，没有对应的 threshold
    valid_idx = np.where(precisions >= (1 - target_fdr))[0]
    
    if len(valid_idx) > 0:
        # 选取第一个满足条件的索引（通常是 PR 曲线中能保证该 Precision 的最低阈值，以最大化 Recall）
        best_idx = valid_idx[0]
        # 边界检查：防止索引超出 thresholds (其长度比 precisions 少 1)
        selected_threshold = thresholds[min(best_idx, len(thresholds)-1)]
    else:
        # 如果没有任何阈值能达到 90% Precision，则退而求其次，选取 Precision 最高的那个
        selected_threshold = thresholds[np.argmax(precisions[:-1])]

    # 3. 基于选定阈值计算具体分类指标
    y_pred = (y_score >= selected_threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

    # 核心指标计算
    tpr = tp / (tp + fn) if (tp + fn) > 0 else 0  # 敏感性 (Sensitivity)
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0  # 误报率 (1 - Specificity)
    actual_precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    actual_fdr = 1 - actual_precision

    return {
        "AUPRC": auprc,
        "AUROC": auroc,
        "TPR_at_10FDR": tpr,           # 核心指标：在误报受控时的检出率
        "FPR_at_10FDR": fpr,           # 核心指标：正常细胞被误诊的概率
        "Actual_FDR": actual_fdr,      # 实际达到的 FDR (应接近 0.1)
        "Actual_Precision": actual_precision, 
        "Threshold": selected_threshold
    }
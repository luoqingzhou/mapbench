from sklearn.metrics import (
    roc_auc_score, 
    average_precision_score, 
    f1_score, 
    precision_recall_curve
)
import numpy as np

def compute_identification_metrics(y_true, y_score, threshold=None):
    """
    支持score分数为任意实数，不限定[0,1]。
    如果threshold为None，将使用F1最优点。
    """
    auroc = roc_auc_score(y_true, y_score)
    auprc = average_precision_score(y_true, y_score)

    precision, recall, thresholds = precision_recall_curve(y_true, y_score)
    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
    best_index = np.argmax(f1_scores)
    best_threshold = thresholds[best_index] if threshold is None else threshold

    y_pred = (y_score >= best_threshold).astype(int)
    f1 = f1_score(y_true, y_pred)

    # return {
    #     "AUROC": auroc,
    #     "AUPRC": auprc,
    #     "F1@best": f1,
    #     "BestThreshold": best_threshold
    # }

    return auprc

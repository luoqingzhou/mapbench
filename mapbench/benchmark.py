from .metrics import metric_template

def evaluate(ref_adata, query_adata, batch_key, label_key, predicted_label_key):
    """
    Run all benchmark metrics.

    Returns:
    - results: dict of metric_name -> score
    """
    results = {}
    results['template_metric'] = metric_template(
        ref_adata, query_adata, batch_key, label_key, predicted_label_key
    )
    return results

def check_anndata(adata, required_keys):
    for key in required_keys:
        if key not in adata.obs:
            raise ValueError(f"Missing key '{key}' in AnnData.obs")

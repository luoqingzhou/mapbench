import os
import numpy as np
import pandas as pd
import anndata as ad
import pytest

from mapbench.metrics import uncertainty

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
COMBINED_PATH = os.path.join(DATA_DIR, "combined_adata.h5ad")

@pytest.fixture
def get_test_data():
    if os.path.exists(COMBINED_PATH):
        print("Using provided combined .h5ad file for testing.")
        adata = ad.read_h5ad(COMBINED_PATH)
    else:
        print("Using synthetic small combined data for testing.")

        # Simulate combined data
        n_ref = 100
        n_query = 20
        n_total = n_ref + n_query
        X = np.random.rand(n_total, 10)

        obs = pd.DataFrame({
            "ref_query": ["ref"] * n_ref + ["query"] * n_query,
            "cell_type": np.random.choice(["A", "B", "C"], size=n_ref).tolist() + [np.nan] * n_query
        }, index=[f"cell_{i}" for i in range(n_total)])

        adata = ad.AnnData(X=X, obs=obs)

    return adata

def test_uncertainty_runs(get_test_data):
    adata = get_test_data

    uncertainties = uncertainty(
        adata=adata,
        embedding="X",
        label_key="cell_type",
        ref_query_key="ref_query",
        ref_key="ref",
        query_key="query",
        n_neighbors=5
    )

    assert isinstance(uncertainties, pd.Series), "uncertainties should be a pandas Series"

    # Check indices match query cells
    query_idx = adata.obs_names[adata.obs["ref_query"] == "query"]
    assert all(uncertainties.index == query_idx), "Indices of uncertainties do not match query cells"
    # Check no NA values
    assert not uncertainties.isna().any(), "uncertainties contains NaN values"

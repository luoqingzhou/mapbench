import numpy as np
from sklearn.neighbors import KNeighborsTransformer
import anndata as ad

# def train_weighted_knn(train_adata: ad.AnnData, embedding: str = "X", n_neighbors: int = 50):
#     """
#     Train a weighted KNN model on reference AnnData.
#     """
#     if embedding == "X":
#         train_emb = train_adata.X
#     elif embedding in train_adata.obsm:
#         train_emb = train_adata.obsm[embedding]
#     else:
#         raise ValueError(f"Embedding '{embedding}' not found in AnnData.")

#     knn_model = KNeighborsTransformer(
#         n_neighbors=n_neighbors,
#         mode="distance",
#         algorithm="brute",
#         metric="euclidean",
#         n_jobs=-1,
#     )
#     knn_model.fit(train_emb)
#     return knn_model

import scanpy as sc
from mapbench import benchmark

# Load example data (replace with your own file paths)
ref_adata = sc.read('reference.h5ad')
query_adata = sc.read('query.h5ad')

results = benchmark.evaluate(
    ref_adata, query_adata,
    batch_key='batch',
    label_key='cell_type',
    predicted_label_key='predicted_cell_type'
)

print("Benchmark results:", results)

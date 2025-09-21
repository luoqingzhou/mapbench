# mapbench

A benchmarking toolkit for evaluating reference mapping in single-cell atlas studies.

## Features
- Evaluate mapping accuracy and anomaly detection
- Modular metrics for biological conservation and abnormal state detection
- Compatible with AnnData (Scanpy) ecosystem

## Installation
1. Configure mapbench Conda Environment

```bash
conda create -n mapbench python=3.10 && y && conda activate mapbench
```

2.  Install Required Environments and Dependencies
Install R environment, compilation tools, and prerequisite packages:

```bash
conda install r-base=4.3.1
conda install compilers
conda install -c conda-forge xz zlib
conda install bioconda::bioconductor-edger=4.0.16
```

3. Compile and Install milopy Package
```bash
cd milopy
pip install .
```

4. Install mapbench Package
```bash
cd ..
pip install .
```
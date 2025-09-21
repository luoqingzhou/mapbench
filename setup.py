from setuptools import setup, find_packages

setup(
    name='mapbench',
    version='0.1.0',
    description='Reference mapping benchmarking toolkit for single-cell atlas',
    author='Zikang Yin',
    author_email='yzk23@mails.tsinghua.edu.cn',
    packages=find_packages(),
     install_requires=[
        'scanpy',
        'anndata',
        'numpy',
        'scikit-learn',
        'pandas',
        'matplotlib',
        'mapqc',
        'jupyter'
    ],
    python_requires='>=3.10',
)

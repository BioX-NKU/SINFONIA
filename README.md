[![PyPI](https://img.shields.io/pypi/v/sinfonia.svg)](https://pypi.org/project/sinfonia)
[![Documentation Status](https://readthedocs.org/projects/sinfonia-svg/badge/?version=latest)](https://sinfonia-svg.readthedocs.io/en/latest/?badge=stable)
[![Downloads](https://pepy.tech/badge/sinfonia)](https://pepy.tech/project/sinfonia)
# SINFONIA: scalable identification of spatially variable genes for deciphering spatial domains


### Find more details on [the Documentation of SINFONIA](https://sinfonia-svg.readthedocs.io/en/latest/index.html).

## Installation
SINFONIA is available on PyPI [here](https://pypi.org/project/sinfonia/) and can be installed via

```
pip install sinfonia
```

You can also install SINFONIA from GitHub via
```
git clone git://github.com/BioX-NKU/SINFONIA.git
cd SINFONIA
python setup.py install
```
The dependencies will be automatically installed along with SINFONIA.


## Quick Start

### Input

* **adata**:       AnnData object of shape `n_obs` × `n_vars`. Rows correspond to cells and columns to genes.
* **mode**:        Mode for identifying spatially variable genes. By default, mode='moran_geary'.
* **n_top_genes**: Number of spatially variable genes to keep. By default, n_top_genes=2000.

### Output

* **adata**:       AnnData object with identified spatially variable genes and additional fields.

SINFONIA can also be seamlessly integrated with [SCANPY](https://scanpy.readthedocs.io/en/stable/), a widely-used Python library for single-cell data analysis.
```python
	import scanpy as sc
	import sinfonia
	# Load the spatial transcriptomic data as an AnnData object (adata)
	# Normalize and logarithmize if the data contains raw counts
	sc.pp.normalize_total(adata)
	sc.pp.log1p(adata)
	# Run SINFONIA
	adata = sinfonia.spatially_variable_genes(adata)
```
### Documentation notebook
We provide a [quick-start notebook](https://github.com/BioX-NKU/SINFONIA/blob/main/docs/source/10X_DLPFC_151507.ipynb) which describes the fundamentals in detail and reproduces the results of SINFONIA. We also provide rich [documentation](https://sinfonia-svg.readthedocs.io/en/latest/index.html) in the form of functional application programming interface documentation, tutorials and example workflows. 

# FlowGrid

FlowGrid density-based clustering algorithm that can perform fast and accurate clustering on very large scRNA-seq data sets. It can be implemented with Scanpy for fast clustering of Scanpy Anndata.

### Installation
FlowGrid supports pip installation.
```sh
pip install sc_FlowGrid
```

### Example1:
Running Flowgrid within Scanpy for scRNA-seq analysis


| requirement | location |
| ------ | ------ |
| Package: Scanpy | https://scanpy.readthedocs.io/en/stable/ |
| Data: Mouse gut endoderm data set | [https://endoderm-explorer.com/] | 
#### Install the packages
```sh
pip install sc_FlowGrid
pip install scanpy
```
#### Import the packages and do the basic setting
```sh
from time import time
import sc_FlowGrid as sf
import scanpy.api as sc
import sys
import warnings
warnings.filterwarnings('ignore')
sc.settings.verbosity = 3  
sc.settings.set_figure_params(dpi=80)  # low dpi (dots per inch) yields small inline figures
```
#### Load the data

```sh
#You can change your file location here
adata = sc.read('~/sc_endoderm_all_cells.h5ad')
```
#### Preprocess
```sh
#Basic filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
#Normalization
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata
#Highly variable genes selection
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var['highly_variable']]
```
#### PCA for dimensionality reduction
```sh
#PCA to 5 dimensions#
sc.tl.pca(adata, n_comps=5)
```
#### Cluster using FlowGrid
You can use sc_autoFlowGrid to do clustering for the data. 
Sample usage is as follows:
```sh
#recomm_parameters = sc_autoFlowGrid(adata, int(set_n), list(binN_range), list(eps_range))
```
Default binN_range = [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25] 
Default eps_range = [1.2, 1.6, 1.9, 2.1, 2.3, 2.7]
By default, 16*6 executions will be performed and the clustering results will be stored in Anndata.obs. By specifing set_n, n recommended parameters will be returned as a list. Users can also specify binN_range and eps_range to reduce computational time such as:
```sh
recomm_parameters = sf.sc_autoFlowGrid(adata, 5, [10, 11, 12, 13, 14, 15, 16, 17, 18, 19], [1.2, 1.6, 1.9, 2.1, 2.3])
```
#### Visualize the result
```sh
#neighbor graph
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=5)
#umap
sc.tl.umap(adata)

#results of recommended parameters
for i in range(len(recomm_parameters)):
    sc.pl.umap(adata, color=recomm_parameters[i],frameon =False)
```

### NOTE
The result of the steps above can be found in the sc_FlowGrid Example.ipynb



### Todos

 - Write MORE Tests
 - Add figures

License
----

MIT


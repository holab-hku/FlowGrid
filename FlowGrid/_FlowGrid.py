import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from typing import Optional
from FlowGrid._FlowGrid_fun import *
bin_n = 9
eps = 1.1
MinDenC = 40  
def FlowGrid(
    adata: AnnData,
    bin_n: int = bin_n,
    eps: float = eps, 
    MinDenC: Optional[int] = None,
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> Optional[AnnData]:
    
    adata = adata.copy() if copy else adata
    projected_data = adata.obsm['X_pca']
    
    n_dimension = projected_data.shape[1]

    transposed_pdata =[]
    for i in range(n_dimension):
        positive_pdata = projected_data[:,i] - np.amin(projected_data[:,i])
        transposed_pdata.append(positive_pdata)
    shifted_pdata = np.array(transposed_pdata).transpose()
    #######################remove dtype############################
    shifted_pdata = shifted_pdata.astype('Float64')    
    if MinDenC:
        fg_object = FlowGrid_fun(shifted_pdata,bin_n = bin_n,eps = eps, MinDenC = MinDenC)
        key_added = 'binN_' + str(bin_n) + '_eps_' + str(eps) + '_MinDenc_' + str(MinDenc) + '_FlowGrid'
    else:
        fg_object = FlowGrid_fun(shifted_pdata,bin_n = bin_n,eps = eps)
        key_added = 'binN_' + str(bin_n) + '_eps_' + str(eps)  + '_FlowGrid'

    flowgrid_label = fg_object.clustering()

    adata.obs[key_added] = pd.Categorical(
    values=flowgrid_label.astype('U'),
    categories=natsorted(np.unique(flowgrid_label).astype('U')),
    )

    return adata if copy else None
            
            

import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from typing import Optional
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.metrics.cluster import v_measure_score
from sklearn.metrics import precision_recall_fscore_support
from FlowGrid._FlowGrid import *
bin_n = 9
eps = 1.1
MinDenB = 10  
MinDenC = 40  
def cluster(
    adata: AnnData,
    MinDenB: Optional[int] = None,
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
    if MinDenB or MinDenC:
        fg_object = FlowGrid(shifted_pdata, MinDenB = MinDenB, bin_n = bin_n, eps = eps, MinDenC = MinDenC)
        key_added = 'MinDenB_'+ str(MinDenB) + '_binN_' + str(bin_n) + '_eps_' + str(eps) + '_MinDenC_' + str(MinDenC) + '_FlowGrid'
    else:
        fg_object = FlowGrid(shifted_pdata, bin_n = bin_n, eps = eps)
        key_added = 'binN_' + str(bin_n) + '_eps_' + str(eps)  + '_FlowGrid'

    flowgrid_label = fg_object.clustering()

    adata.obs[key_added] = pd.Categorical(
    values=flowgrid_label.astype('U'),
    categories=natsorted(np.unique(flowgrid_label).astype('U')),
    )

    return adata if copy else None

labels_pred = ['louvain']
labels_ref = ['Cluster']
def AdjustedRandScore(
    adata: AnnData,
    labels_pred: list = labels_pred,
    labels_ref: list = labels_ref, 
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> Optional[AnnData]:
    
    adata = adata.copy() if copy else adata

    df = adata.obs[labels_pred + labels_ref]
    
    for labelP in labels_pred:
        for labelR in labels_ref:
            df_filtered  = df[df[labelP] != str(-1.0)]
            df_filtered = df_filtered.applymap(str)
            print(labelP +" vs." + labelR + " ARI:"+ str(round(ARI(df_filtered[labelP].values, df_filtered[labelR].values),4)))

remain_list = ['binN_12_eps_1.2_FlowGrid']
def keep_labels(
    adata: AnnData,
    remain_list: list = remain_list,
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> Optional[AnnData]:
    adata = adata.copy() if copy else adata
    labels = adata.obs.columns
    for label in labels:
        if 'binN_' in label:
            if label in remain_list:
                continue
            else:
                del adata.obs[label]
        else:
            continue          

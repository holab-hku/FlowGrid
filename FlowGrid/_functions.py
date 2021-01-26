import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from typing import Optional
from sklearn.metrics import adjusted_rand_score as ARI
from sklearn.metrics.cluster import v_measure_score
from sklearn.metrics import precision_recall_fscore_support
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

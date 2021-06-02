import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from typing import Optional
from sklearn.metrics import adjusted_rand_score
from FlowGrid._FlowGrid import *


bin_n = 15
eps = 5
MinDenB = 200
MinDenC = 500 
def cluster(
    adata: AnnData,
    MinDenB: Optional[int] = 200,
    bin_n: int = bin_n,
    eps: float = eps, 
    MinDenC: Optional[int] = 500,
    DimRMethod: Optional[str] = 'PCA',
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,    
) -> Optional[AnnData]:
    
    adata = adata.copy() if copy else adata
    if DimRMethod == 'PCA':
        projected_data = adata.obsm['X_pca'][:,0:5].astype(float)
    if DimRMethod == 'UMAP':
        projected_data = adata.obsm['X_umap'][:,0:5].astype(float)
    if DimRMethod == 'DIFFMAP':
        projected_data = adata.obsm["X_diffmap"][:,0:5].astype(float)
    if DimRMethod == 'NMF':
        projected_data = adata.obsm['X_nmf'][:,0:5].astype(float)
    if DimRMethod == 'FA':
        projected_data = adata.obsm['X_fa'][:,0:5].astype(float)   

        
    projected_data.shape[1]
    n_dimension = 5
    transposed_pdata =[]
    for i in range(n_dimension):
        positive_pdata = projected_data[:,i] - np.amin(projected_data[:,i])
        transposed_pdata.append(positive_pdata)
    shifted_pdata = np.array(transposed_pdata).transpose()
    shifted_pdata = shifted_pdata.astype(float)

    fg_object = FlowGrid(shifted_pdata, MinDenB = MinDenB, bin_n = bin_n, eps = math.sqrt(eps), MinDenC = MinDenC)
    key_added = 'MinDenB_'+ str(MinDenB) + '_binN_' + str(bin_n) + '_eps_' + str(eps) + '_MinDenC_' + str(MinDenC) + '_FlowGrid'

    flowgrid_label = fg_object.clustering()

    adata.obs[key_added] = pd.Categorical(
    values=flowgrid_label.astype('U'),
    categories=natsorted(np.unique(flowgrid_label).astype('U')),
    )

    return adata if copy else None
            
remain_list = ['FlowGrid1']
def clean_labels(
    adata: AnnData,
    remain_list: list = remain_list,
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> Optional[AnnData]:
    adata = adata.copy() if copy else adata
    labels = adata.obs.columns
    flabels = []
    for i in range(len(remain_list)):
        flabels.append('FlowGrid' + str(i+1))
    adata.obs[flabels] = adata.obs[remain_list]
    for label in adata.obs.columns:
        if 'binN_' in label:
            del adata.obs[label]
            
labels_pred = ['louvain']
labels_ref = ['Cluster']
def AdjustedRandScore(
    adata: AnnData,
    labels_pred: list = labels_pred,
    labels_ref: list = labels_ref, 
    removenoise: bool = False,
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> Optional[AnnData]:
    adata = adata.copy() if copy else adata
    df = adata.obs[labels_pred + labels_ref]
    ARI_values = {}
    for labelP in labels_pred:
        ARI_values[labelP] = []
        for labelR in labels_ref:            
            if removenoise == True:
                df_filtered  = df[df[labelP] != str(-1.0)]
                df_filtered  = df_filtered[df_filtered[labelR] != 'nan']
                df_filtered = df_filtered.applymap(str)
                
                ARI_values[labelP].append(round(adjusted_rand_score(df_filtered[labelP].values, df_filtered[labelR].values),2))
            else:
                df = df.applymap(str)
                ARI_values[labelP].append(round(adjusted_rand_score(df[labelP].values, df[labelR].values),2))
    ARIinfo = pd.DataFrame.from_dict(ARI_values).T
    ARIinfo.columns = labels_ref
    return ARIinfo        

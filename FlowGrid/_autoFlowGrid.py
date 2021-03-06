import numpy as np
import math
from sklearn.preprocessing import LabelEncoder
from anndata import AnnData
from typing import Optional
from time import time
from FlowGrid._functions import *

def calinski_harabasz_score(X, labels):
    """Compute the Calinski and Harabasz score.
    It is also known as the Variance Ratio Criterion.
    The score is defined as ratio between the within-cluster dispersion and
    the between-cluster dispersion.
    Read more in the :ref:`User Guide <calinski_harabasz_index>`.
    Parameters
    ----------
    X : array-like, shape (``n_samples``, ``n_features``)
        List of ``n_features``-dimensional data points. Each row corresponds
        to a single data point.
    labels : array-like, shape (``n_samples``,)
        Predicted labels for each sample.
    Returns
    -------
    score : float
        The resulting Calinski-Harabasz score.
    References
    ----------
    .. [1] `T. Calinski and J. Harabasz, 1974. "A dendrite method for cluster
       analysis". Communications in Statistics
       <https://www.tandfonline.com/doi/abs/10.1080/03610927408827101>`_
    """
    le = LabelEncoder()
    labels = le.fit_transform(labels)

    n_samples, _ = X.shape
    n_labels = len(le.classes_)

    extra_disp, intra_disp = 0., 0.
    mean = np.mean(X, axis=0)
    for k in range(n_labels):
        cluster_k = X[labels == k]
        mean_k = np.mean(cluster_k, axis=0)
        extra_disp += len(cluster_k) * np.sum((mean_k - mean) ** 2)
        intra_disp += np.sum((cluster_k - mean_k) ** 2)

    return (1. if intra_disp == 0. else
            extra_disp * (n_samples - n_labels) /
            (intra_disp * (n_labels - 1.)))

def autoFlowGrid(
    adata: AnnData,
    nTop: Optional[int] = 5,
    MinDenB: Optional[int] = 200,
    Bin_n: Optional[list] = None,    
    Eps: Optional[list] = None,
    MinDenC: Optional[int] = 500,
    nCluster: Optional[int] = 5,
    penalty: Optional[float] = 1.0,
    MinCHI: Optional[int] = 5000,
    MaxNoise: Optional[float] = 0.15,
    DimRMethod: Optional[str] = 'PCA',
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None
    
) -> Optional[AnnData]:
    t0=time()
    adata = adata.copy() if copy else adata

    ###################################################################
    Bin_n = Bin_n if Bin_n else [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
    Eps = Eps if Eps else [1,2,3,4,5]
    if DimRMethod == 'PCA':
        feature_data = adata.obsm['X_pca'].astype(float)
    if DimRMethod == 'UMAP':
        feature_data = adata.obsm['X_umap'].astype(float)
    if DimRMethod == 'DIFFMAP':
        feature_data = adata.obsm["X_diffmap"].astype(float)
    if DimRMethod == 'NMF':
        feature_data = adata.obsm['X_nmf'].astype(float)
    if DimRMethod == 'FA':
        feature_data = adata.obsm['X_fa'].astype(float)    
        
    if len(adata) < 100000:
        MinDenB = 50
        MinDenC = 100
        
    CHixNobs_values = {}
    
    for eps in Eps:
        MaxNobs = -1       
        for bin_n in Bin_n:
            
            cluster(adata,MinDenB,bin_n,eps,MinDenC,DimRMethod)
                                  
            label_data = adata.obs['MinDenB_'+ str(MinDenB)+'_binN_'+str(bin_n)+'_eps_'+ str(eps)+ '_MinDenC_' + str(MinDenC)+'_FlowGrid'].tolist()
            if MaxNobs == -1:
                MaxNobs = len(set(label_data))
            else:
                if len(set(label_data)) >= MaxNobs:
                    MaxNobs = len(set(label_data))
                else:
                    MaxNobs = MaxNobs
            if label_data.count('-1.0')/len(label_data) > MaxNoise:
                break
            if len(set(label_data)) < nCluster:
                del adata.obs['MinDenB_'+ str(MinDenB)+'_binN_'+str(bin_n)+'_eps_'+str(eps)+'_MinDenC_' + str(MinDenC)+'_FlowGrid']
                continue
            indices = [i for i, x in enumerate(label_data) if x != "-1.0"]
            CHI = calinski_harabasz_score(feature_data[indices], list(filter(lambda a: a != '-1.0', label_data)))
            if CHI < MinCHI or CHI != CHI:
                del adata.obs['MinDenB_'+ str(MinDenB)+'_binN_'+str(bin_n)+'_eps_'+str(eps)+'_MinDenC_' + str(MinDenC)+'_FlowGrid']
                continue

            CHixNobs_values['MinDenB_'+ str(MinDenB)+'_binN_'+str(bin_n)+'_eps_'+str(eps)+'_MinDenC_' + str(MinDenC)+'_FlowGrid']\
            = [round(CHI * len(set(label_data))**(penalty),0), round(CHI,0), len(set(label_data)), label_data.count('-1.0')/len(label_data)] 
        if MaxNobs < nCluster:
            break            
            
    maxn_parameters = sorted(CHixNobs_values, key=CHixNobs_values.get, reverse=True)[:nTop]
    CHIinfo = pd.DataFrame.from_dict(CHixNobs_values).T
    CHIinfo = CHIinfo.rename(columns={0: 'uCHI',1: 'CHI',2: 'Nobs',3: 'OurliersPortion'})
    
    if nTop <= len(CHixNobs_values):
        print("autoFlowGrid completed in : "+ str(round(time()-t0,3)) + " seconds.\n"\
          + str(len(CHixNobs_values)) + " sets of parameters are stored.\n"\
          + str(nTop) + " sets of parameters are recommended.\n")
    else:
        print("autoFlowGrid completed in : "+ str(round(time()-t0,3)) + " seconds.\n"\
          + "There are only " + str(len(CHixNobs_values)) + " potential good parameters.\n"\
          +  str(len(CHixNobs_values)) + " sets of parameters are recommended.\n")
    return maxn_parameters, CHIinfo
    
    
    

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
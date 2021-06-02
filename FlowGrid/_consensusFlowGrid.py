import ClusterEnsembles as CE 
import numpy as np
import pandas as pd
from anndata import AnnData
from natsort import natsorted
from typing import Optional
from time import time
from FlowGrid._FlowGrid import *

def consensusFlowGrid(
    adata: AnnData,
    nTop: Optional[int] = 5,
    MinDenB: Optional[int] = None,
    Bin_n: Optional[list] = None,    
    Eps: Optional[list] = None,
    MinDenC: Optional[int] = None,
    nCluster: Optional[int] = 5,
    penalty: Optional[float] = 1.0,
    MinCHI: Optional[int] = 5000,
    MaxNoise: Optional[float] = 0.15,
    DimRMethod: Optional[str] = 'PCA',
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
    nDims: Optional[int] = 20,
    nBases: Optional[int] = None,
    nClass: Optional[int] = None,
    
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
        
    feature_data2 = feature_data[:,(0,1,2,3,4)]

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
            
    maxn_parameters1 = sorted(CHixNobs_values, key=CHixNobs_values.get, reverse=True)[:1]
    ps = maxn_parameters1[0].split('_')
    adata.obs['autoFlowGrid'] = adata.obs[maxn_parameters1]
    dims= adata.obs[['autoFlowGrid']]
    dims.columns = [0]
    if DimRMethod == 'PCA':
        projected_data = adata.obsm['X_pca'].astype(float)
    if DimRMethod == 'UMAP':
        projected_data = adata.obsm['X_umap'].astype(float)
    if DimRMethod == 'DIFFMAP':
        projected_data = adata.obsm["X_diffmap"].astype(float)
    if DimRMethod == 'NMF':
        projected_data = adata.obsm['X_nmf'].astype(float)
    if DimRMethod == 'FA':
        projected_data = adata.obsm['X_fa'].astype(float) 
    list_dims = list(range(3,nDims))
    np.random.shuffle(list_dims)
    if nBases:
        for r in range(nBases):
            rdims3 = random.sample(range(3, nDims), 3)
            projected6_data = projected_data[:,tuple(rdims3 + [0,1,2])]    
            CHixNobs_values = {}
            n_dimension = projected5_data.shape[1]
            transposed_pdata =[]
            for i in range(n_dimension):
                positive_pdata = projected6_data[:,i] - np.amin(projected6_data[:,i])
                transposed_pdata.append(positive_pdata)
            shifted_pdata = np.array(transposed_pdata).transpose()
            shifted_pdata = shifted_pdata.astype(float)
            fg_object = FlowGrid(shifted_pdata, MinDenB = int(ps[1]), bin_n = int(ps[3]), eps = math.sqrt(float(ps[5])), MinDenC = int(ps[7]))
            flowgrid_label = fg_object.clustering()
            dims[r + 1] = flowgrid_label
    else:
        for r in range((len(range(3,nDims))) // 3):
            dims3, list_dims = list_dims[-3:], list_dims[0:-3]

            projected6_data = projected_data[:,tuple(dims3 + [0,1,2])]    
            CHixNobs_values = {}
            n_dimension = projected6_data.shape[1]
            transposed_pdata =[]
            for i in range(n_dimension):
                positive_pdata = projected6_data[:,i] - np.amin(projected6_data[:,i])
                transposed_pdata.append(positive_pdata)
            shifted_pdata = np.array(transposed_pdata).transpose()
            shifted_pdata = shifted_pdata.astype(float)
            fg_object = FlowGrid(shifted_pdata, MinDenB = int(ps[1]), bin_n = int(ps[3]), eps = math.sqrt(float(ps[5])), MinDenC = int(ps[7]))
            flowgrid_label = fg_object.clustering()
            dims[r + 1] = flowgrid_label
    if nClass:
        consF = CE.cluster_ensembles(dims.values.T, nclass=nClass,solver='mcla',verbose = True)
    else:  
        consF = CE.cluster_ensembles(dims.values.T, nclass=len(set(adata.obs['autoFlowGrid'])),solver='mcla',verbose = True)
    adata.obs['consensusFlowGrid'] = consF.astype(str)
    for label in adata.obs.columns:
        if 'binN_' in label:
            del adata.obs[label]
    print("ConsensusFlowGrid completed in : "+ str(round(time()-t0,3)) + " seconds.\n")
    
    

import numpy as np
from sklearn.preprocessing import LabelEncoder
from anndata import AnnData
from typing import Optional
from time import time
#from FlowGrid._functions import *
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
    #X, labels = check_X_y(X, labels)
    le = LabelEncoder()
    labels = le.fit_transform(labels)

    n_samples, _ = X.shape
    n_labels = len(le.classes_)

    #check_number_of_labels(n_labels, n_samples)

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

set_n = 5
MinDenB = 10
MinDenC = 40
penalty = 1.0
def autoFlowGrid(
    adata: AnnData,
    set_n: int = set_n,
    MinDenB: Optional[int] = None,
    Bin_n: Optional[list] = None,    
    Eps: Optional[list] = None,
    MinDenC: Optional[int] = None,
    penalty: Optional[int] = 1.0,
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None
) -> Optional[AnnData]:
    t0=time()
    print("autoFlowGrid Starts!")
    adata = adata.copy() if copy else adata
    ###################################################################
    Bin_n = Bin_n if Bin_n else [5,6,7,8,9,10,11,12,13,14,15]
    Eps = Eps if Eps else [0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
    feature_data = adata.obsm['X_pca'].astype('Float64')
    CHixNobs_values = {}
    
    for eps in Eps:
        maxbin_n = 24
        for bin_n in Bin_n:            
            cluster(adata,MinDenB,bin_n,eps,MinDenC)
            if MinDenB or MinDenC:
                label_data = adata.obs['MinDenB_'+ str(MinDenB)+'_binN_'+str(bin_n)+'_eps_'+ str(eps)+ '_MinDenC_' + str(MinDenC)+'_FlowGrid'].tolist()
                if label_data.count('-1.0')/len(label_data) > 0.2:
                    maxbin_n = bin_n
                    break
                CHI = calinski_harabasz_score(feature_data, label_data)
                if CHI < 10000 or CHI != CHI:
                    del adata.obs['MinDenB_'+ str(MinDenB)+'_binN_'+str(bin_n)+'_eps_'+str(eps)+'_MinDenC_' + str(MinDenC)+'_FlowGrid']
                    continue
                if len(set(label_data)) < 6 or len(set(label_data)) > 60:
                    del adata.obs['MinDenB_'+ str(MinDenB)+'_binN_'+str(bin_n)+'_eps_'+str(eps)+'_MinDenC_' + str(MinDenC)+'_FlowGrid']
                    continue
                CHixNobs_values['MinDenB_'+ str(MinDenB)+'_binN_'+str(bin_n)+'_eps_'+str(eps)+'_MinDenC_' + str(MinDenC)+'_FlowGrid']\
                = [round(CHI * len(set(label_data))**(penalty)), round(CHI), int(len(set(label_data)))] 
            else:
                label_data = adata.obs['binN_'+str(bin_n)+'_eps_'+ str(eps)+'_FlowGrid'].tolist()
                if label_data.count('-1.0')/len(label_data) > 0.2:
                    maxbin_n = bin_n
                    break
                CHI = calinski_harabasz_score(feature_data, label_data)
                if CHI < 10000 or CHI != CHI:
                    del adata.obs['binN_'+str(bin_n)+'_eps_'+str(eps)+'_FlowGrid']
                    continue
                if len(set(label_data)) < 6 or len(set(label_data)) > 60:
                    del adata.obs['binN_'+str(bin_n)+'_eps_'+str(eps)+'_FlowGrid']
                    continue                                        
                CHixNobs_values['binN_'+str(bin_n)+'_eps_'+str(eps)+'_FlowGrid'] \
                = [round(CHI * len(set(label_data))**(penalty)),round(CHI),int(len(set(label_data)))]

    CHI_report = pd.DataFrame.from_dict(CHixNobs_values).T.rename(columns={0: "iCHI", 1: "CHI", 2:"nobs"})
    a = CHI_report.drop_duplicates(subset =["iCHI"], keep = 'first', inplace = False)   
    recom_parameters = a.sort_values(by=['iCHI'],ascending=False).index[0:set_n].tolist()
    CHI_report = a.sort_values(by=['iCHI'],ascending=False)
    #maxn_CHixNobs_parameters = sorted(CHixNobs_values, key=CHixNobs_values.get, reverse=True)[:set_n]

    if set_n < len(CHixNobs_values):
        print("autoFlowGrid completed in : "+ str(round(time()-t0,3)) + " seconds.\n"\
          + str(len(CHixNobs_values)) + " sets of parameters are stored.\n"\
          + str(set_n) + " sets of parameters are recommended.\n")
    else:
        print("autoFlowGrid completed in : "+ str(round(time()-t0,3)) + " seconds.\n"\
          + "There are only " + str(len(CHixNobs_values)) + " potential good results.\n"\
          +  str(len(CHixNobs_values)) + " sets of parameters are recommended.\n")
    return recom_parameters, CHI_report
    
    
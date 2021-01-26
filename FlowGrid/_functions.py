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
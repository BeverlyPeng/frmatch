
import numpy as np
import pandas as pd
from sklearn.preprocessing import normalize
from scipy.sparse import csr_matrix
from tqdm import tqdm
import time

def normalization(adata, cluster_header, scale = True, norm_by = None, save = False, verbose = 0): 
    """\
    Normalizing cell by gene data in `adata` using pd.DataFrame.

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        scale: bool (default: True)
            Whether to scale to 0 and 1 per gene.
        norm_by: ["median", "mean"] (default: None)
            How to normalize cell by gene.
        save: bool | str (default: False)
            Whether to save results as csv. If string, save as string.
        verbose: int (default: 0)
            Controls print statements.
    
    Returns
    -------
    adata: AnnData
        AnnData with normalized data stored in adata.X
    """
    start_time = time.time()
    if verbose > 0: print("--- %s seconds ---" % (time.time() - start_time))

    mat = adata.to_df().copy() # pandas
    
    ## scale to 0 and 1 PER GENE
    if scale: 
        # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.div.html
        mat = mat.div(mat.max(axis=0), axis=1) #### col-wise normalizing
    # concurrent.futures.ProcessPoolExecution
    if verbose > 0: print("--- %s seconds ---" % (time.time() - start_time))
    clusters = list(np.unique(adata.obs[cluster_header]))
    for cluster in tqdm(np.unique(clusters), desc = "normalizing"): 
        curr_time = time.time()
        subset = mat.loc[adata.obs[cluster_header] == cluster,:].copy()
        mat = mat.drop(subset.index, axis = 0)
        ## normalize by row median
        if norm_by == "median": 
            subset = subset.mul(subset.median(axis=0), axis=1) 
        ## normalize by row mean
        elif norm_by == "mean": 
            subset = subset.mul(subset.mean(axis=0), axis=1) 
        ## final re-scaling to 0 and 1
        subset = subset / subset.max().max()
        mat = pd.concat([mat, subset])
        if save == True: 
            mat.to_csv("normalized_matrix_df.csv", index = True)
        elif isinstance(save, str): 
            mat.to_csv(save, index = True)
        del subset
        if verbose > 0: print("--- %s seconds ---" % (time.time() - curr_time), cluster)
        
    # Sorting to original sample order
    mat = mat.iloc[pd.Categorical(mat.index, adata.obs.index).argsort()].replace(float('nan'), 0.0)
    adata.X = mat
    if verbose > 0: print("--- %s seconds ---" % (time.time() - start_time))
    return adata

def normalization_np(adata, cluster_header, scale = True, norm_by = None, save = None, verbose = 0): 
    """\
    Normalizing cell by gene data in `adata` using numpy.

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        scale: bool (default: True)
            Whether to scale to 0 and 1 per gene.
        norm_by: ["median", "mean"] (default: None)
            How to normalize cell by gene.
        save: bool | str (default: False)
            Whether to save results as csv. If string, save as string.
        verbose: int (default: 0)
            Controls print statements.
    
    Returns
    -------
    adata: AnnData
        AnnData with normalized data stored in adata.X
    """
    # data = np.array([[1,1,1],[2,2,2],[3,3,3],[4,4,4]]) # (4, 3)
    # # data[0] # [1, 1, 1]
    # data = data / np.amax(data, axis = 0) # max per col [4, 4, 4]
    # data = data / np.amax(data, axis = 1)[:, None] # max per row [1, 2, 3, 4]
    # data = data * np.median(data, axis = 0) # [0.625, 0.625, 0.625]
    start_time = time.time()
    if verbose > 0: print("--- %s seconds ---" % (time.time() - start_time))
    mat = adata.to_df().copy().to_numpy() # numpy
    # ## scale to 0 and 1 PER GENE
    if scale: 
        # https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.div.html
    #     mat = mat / np.amax(mat, axis = 0)
        mat = np.divide(mat, np.amax(mat, axis = 0), out = np.zeros_like(mat), where=np.amax(mat, axis = 0)!=0)
    # concurrent.futures.ProcessPoolExecution
    if verbose > 0: print("--- %s seconds ---" % (time.time() - start_time))
    clusters = list(np.unique(adata.obs[cluster_header]))
    mat_norm = np.empty(mat.shape)
    # cluster = "e1_e299_SLC17A7_L5b_Cdh13"

    for cluster in tqdm(np.unique(clusters), desc = "normalizing"): 
    #     print(cluster)
        curr_time = time.time()
        keep = [i for i in range(len(list(adata.obs[cluster_header]))) if list(adata.obs[cluster_header])[i] == cluster]
    #     print(cluster, len(keep), keep[:10])
        subset = mat[keep]
        ## normalize by row median
        if norm_by == "median": 
            subset = subset * np.median(subset, axis = 0)
        ## normalize by row mean
        elif norm_by == "mean": 
            subset = subset * np.mean(subset, axis = 0)

        ## final re-scaling to 0 and 1
        subset = subset / max(np.amax(subset, axis = 0))
    #     subset = np.divide(subset, max(np.amax(subset, axis = 0)), out = np.zeros_like(subset), where=np.amax(subset, axis = 0)!=0)
        mat_norm[keep] = subset
        if verbose > 0: print("--- %s seconds ---" % (time.time() - curr_time), cluster)
        
    # Sorting to original sample order
    # mat = mat.iloc[pd.Categorical(mat.index, adata.obs.index).argsort()]
    if save == True: 
        np.to_csv("normalized_matrix_np.csv", mat_norm, delimiter = ",")
    elif isinstance(save, str): 
        np.to_csv(save, mat_norm, delimiter = ",")
    adata.X = mat_norm
    if verbose > 0: print("--- %s seconds ---" % (time.time() - start_time))
    return adata

def csr_vappend(a,b):
    """\
    Stack arrays horizontally (column-wise) using numpy.hstack. 
    """
    a.data = np.hstack((a.data,b.data))
    a.indices = np.hstack((a.indices,b.indices))
    a.indptr = np.hstack((a.indptr,(b.indptr + a.nnz)[1:]))
    a._shape = (a.shape[0]+b.shape[0],b.shape[1])
    return a

def normalization_csr(adata, cluster_header, scale = True, norm_by = None, save = False, verbose = 0): 
    """\
    Normalizing cell by gene data in `adata` using csr_matrix.

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        scale: bool (default: True)
            Whether to scale to 0 and 1 per gene.
        norm_by: ["median", "mean"] (default: None)
            How to normalize cell by gene.
        save: bool | str (default: False)
            Whether to save results as csv. If string, save as string.
        verbose: int (default: 0)
            Controls print statements.
    
    Returns
    -------
    adata: AnnData
        AnnData with normalized data stored in adata.X
    """
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
    mat = csr_matrix(adata.to_df().to_numpy())
    ## scale to 0 and 1 PER GENE (PER COL)
    if scale: 
        mat = normalize(mat, axis=0, norm='max') #### col-wise normalizing ## need to double check
    
    clusters = list(np.unique(adata.obs[cluster_header]))
    cluster_assignments = list(adata.obs[cluster_header])
    first = True
    csr_comb = None
    for cluster in clusters: 
        print(cluster)
        subset = [i for i, x in enumerate(cluster_assignments) if x == cluster]
        subset = mat[subset]
        # todo: implement median
        mean = [float(val) for val in subset.mean(axis = 1)]
        ct = 0
        for i in subset.nonzero()[0]: 
            # print(i, subset.data[ct], "x", mean[i], "=", subset.data[ct] / mean[i])
            subset.data[ct] = subset.data[ct] / mean[i]
            ct += 1
        ## final re-scaling to 0 and 1
        subset = subset / subset.max()
        if first: 
            csr_comb = subset.copy()
            first = False
        else: 
            csr_comb = csr_vappend(csr_comb, subset)
        if save == True: 
            mat.to_csv("normalized_matrix_csr.csv", index = True)
        elif isinstance(save, str): 
            mat.to_csv(save, index = True)
        # del subset
        
    # Sorting to original sample order
    # mat = mat.iloc[pd.Categorical(mat.index, adata.obs.index).argsort()]
    # adata.X = mat
    return csr_comb

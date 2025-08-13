
from IPython.core.debugger import set_trace
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy import stats
from scipy.stats import norm
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from statsmodels.stats.multitest import multipletests

def filter_cluster(adata, cluster_header, filter_size): 
    """\
    Filtering out clusters smaller than `filter_size`. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in `adata.obs` storing cell annotation.
        filter_size: float (default: 0.99)
            Minimum cluster size.
    
    Returns
    -------
    adata: AnnData
        Subset AnnData based on `filter_size`.
    """
    var_names = list(adata.var.index)
    # Getting cluster sizes
    tab = pd.DataFrame(adata.obs[cluster_header].value_counts()).reset_index()
    # Filtering out small clusters
    columns = list(tab.columns)
    cluster_keep = list(tab[tab[columns[1]] > filter_size][columns[0]])
    adata = adata[adata.obs[cluster_header].isin(cluster_keep)]
    adata.var.index = var_names
    return adata

def padj_FRmatch(pmat, p_adj_method = "fdr_by"): 
    """\
    P-value adjustment statsmodels.stats.multitest.multipletests.

    Parameters
    ----------
        pmat: pd.DataFrame
            P-value matrix.
        p_adj_method: str (default: "fdr_by")
            P-value adjustment method. 
    
    Returns
    -------
    pvals_adj: pd.DataFrame
        P-value adjusted matrix.
    """
    # https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
    # statsmodels.stats.multitest
    pvals_adj = multipletests(pmat.unstack(), method = p_adj_method)
    pvals_adj = pvals_adj[1].reshape((pmat.shape[1], pmat.shape[0])).T
    pvals_adj = pd.DataFrame(pvals_adj)
    pvals_adj.columns = pmat.columns
    pvals_adj.index = pmat.index
    return pvals_adj

def cutoff_FRmatch(pmat, p_adj_method = "fdr_by", sig_level = 0.05): 
    """\
    Determines matches by cutting off the p-values from FRmatch and adds the "unassigned" row at the bottom.

    Parameters
    ----------
        pmat: pd.DataFrame
            P-value matrix.
        p_adj_method: str (default: "fdr_by")
            P-value adjustment method. 
        sig_level: float (default: 0.05)
            P-value cutoff threshold
    
    Returns
    -------
    pvals_adj: pd.DataFrame
        Matrix where 1 means p-value exceeded `sig_level`.
    """
    # https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
    # pvals_adj = multipletests(pmat.unstack(), method = p_adj_method)
    # pvals_adj = pvals_adj[1].reshape((pmat.shape[1], pmat.shape[0])).T
    # pvals_adj = pvals_adj > sig_level
    # pvals_adj = pd.DataFrame(pvals_adj)
    # pvals_adj.columns = pmat.columns
    # pvals_adj.index = pmat.index

    pvals_adj = padj_FRmatch(pmat, p_adj_method)
    pvals_adj = pvals_adj > sig_level
    pvals_adj = pvals_adj.astype(int)
    zeros = []
    for col in pvals_adj.columns: 
        if sum(pvals_adj[col]) == 0: 
            zeros.append(col)
    temp = pd.DataFrame(dict(zip(pvals_adj.columns, [1 if val in zeros else 0 for val in pvals_adj.columns])), index = ["unassigned"])
    if isinstance(temp.columns, pd.MultiIndex): 
        temp.columns.names = ["query_cluster", "index"]
    else: 
    # if len(list(temp.columns)[0]) == 1: 
        temp.columns.name = "query_cluster"
    # elif len(list(temp.columns)[0]) == 2: 
    #     temp.columns.names = ["query_cluster", "index"]
    # else: 
    #     print("Something wrong")
    #     print(list(temp.columns))
    #     print(temp.columns)
    #     print(temp)
    
        # return

    pvals_adj = pd.concat([pvals_adj, temp])
    return pvals_adj

def reorder_FRmatch(df, axis = 1): 
    """\
    Reorders cutoff_FRmatch output to form a diagonal.

    Parameters
    ----------
        df: pd.DataFrame
            Output of cutoff_FRmatch
        axis: [0, 1] (default: 1)
            Axis to reorder pd.DataFrame.
    
    Returns
    -------
    df: pd.DataFrame
        Ordered pd.DataFrame with a diagonal.
    """
    if axis == 0: 
        df = df.sort_values(by = list(df.columns), axis = 0, ascending = False)
    elif axis == 1: 
        df = df.sort_values(by = list(df.index), axis = 1, ascending = False)
    else: print("Error: `axis` must be in [0, 1]")
    
    return df

def get_annotation(results, p_adj_method = "fdr_by", sig_level = 0.05, force_match = False):
    """Converts frmatch.FRmatch_cell2cluster output into annotation format for adata.obs.
    Applied p-value adjustment and p-value threshold cutoffs.
    Adds "unassigned" for cells with no p-values above cutoff.
    This method is adapted from frmatch.plot_FRmatch_cell2cluster.py.

    Args:
        results (pd.DataFrame): frmatch.FRmatch_cell2cluster output

    Returns:
        pd.DataFrame: Annotation dataframe to be merged with adata.obs.
    """

    # Pivotting results to ref_cluster as index and query_cluster as columns
    results = results.sort_values(
        ["query_cluster", "index", "p_value"], ascending=[True, True, False]
    )
    results = results[
        ~results.duplicated(["index", "query_cluster", "ref_cluster"], keep="first")
    ]
    results = results.pivot(
        index="ref_cluster", columns=["query_cluster", "index"], values="p_value"
    ).replace(np.nan, 0)

    if force_match: 
        max_idx = results.idxmax()
        for col in results.columns: 
            row = max_idx[col]
            results.loc[results.index != row, col] = 0

    # Adjusting p-values
    results_2 = padj_FRmatch(results, p_adj_method=p_adj_method)
    
    # Adding the unassigned row
    if not force_match: 
        # Setting values to zero unless passes sig_level threshold
        results_3 = results_2.applymap(lambda x: x if x > sig_level else 0)

        # Adding the unassigned row
        samples = results_3.columns[results_3.sum() == 0]
        temp = pd.DataFrame(
            dict(zip(samples, [1] * len(samples))), index=["unassigned"]
        )
        if temp.shape[1] != 0:
            temp.columns.names = ["query_cluster", "index"]
            temp.index.name = "ref_cluster"
            results_4 = pd.concat([results_3, temp])
        else:
            results_4 = results_3.copy()
    else: results_4 = results_2.copy()

    # idxmax: returns ref_cluster with largest value per sample
    temp = pd.DataFrame(results_4.max())
    annotation = pd.DataFrame(results_4.idxmax()).reset_index()
    annotation[1] = list(temp[0])
    annotation = annotation.rename(columns = {0: "ref_cluster", 1: "p_value"})

    return annotation

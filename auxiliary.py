
from IPython.core.debugger import set_trace
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx as nx
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib

from statsmodels.stats.multitest import multipletests
import seaborn as sns
from scipy import stats

def filter_cluster(adata, cluster_header, filter_size): 
    """\
    Filtering out clusters smaller than `filter_size`. 
    """
    # Getting cluster sizes
    tab = pd.DataFrame(adata.obs[cluster_header].value_counts()).reset_index()
    # Filtering out small clusters
    columns = list(tab.columns)
    cluster_keep = list(tab[tab[columns[1]] > filter_size][columns[0]])
#     cluster_keep = list(tab[tab["count"] > filter_size][cluster_header])
#     cluster_keep = list(tab[tab[cluster_header] > filter_size]["index"])
    adata = adata[adata.obs[cluster_header].isin(cluster_keep)]
    
    return adata

def padj_FRmatch(pmat, p_adj_method = "fdr_by"): 
    # https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
    # statsmodels.stats.multitest
    pvals_adj = multipletests(pmat.unstack(), method = p_adj_method)
    pvals_adj = pvals_adj[1].reshape((pmat.shape[1], pmat.shape[0])).T
    
#     pvals_adj = stats.false_discovery_control(pmat, method = p_adj_method)
#     pvals_adj = pvals_adj.reshape((pmat.shape[0], pmat.shape[1]))
    
    pvals_adj = pd.DataFrame(pvals_adj)
    pvals_adj.columns = pmat.columns
    pvals_adj.index = pmat.index
    return pvals_adj

def cutoff_FRmatch(pmat, p_adj_method = "fdr_by", sig_level = 0.05): 
    # https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
    pvals_adj = multipletests(pmat.unstack(), method = p_adj_method)
    pvals_adj = pvals_adj[1].reshape((pmat.shape[1], pmat.shape[0])).T
    pvals_adj = pvals_adj > sig_level
    
    # https://scipy.github.io/devdocs/reference/generated/scipy.stats.false_discovery_control.html
#     pvals_adj = stats.false_discovery_control(pmat, method = p_adj_method)
#     pvals_adj = pvals_adj > sig_level
#     pvals_adj = pvals_adj.reshape((pmat.shape[0], pmat.shape[1]))

    pvals_adj = pd.DataFrame(pvals_adj)
    pvals_adj.columns = pmat.columns
    pvals_adj.index = pmat.index
    print(pvals_adj.columns)
    print(pvals_adj.index)
    pvals_adj = pvals_adj.astype(int)
    pvals_adj = pvals_adj.astype(int)
    # pvals_adj = pd.concat([pvals_adj, pd.DataFrame(dict(pvals_adj.sum(axis = 0)), index = ["unassigned"])])
    zeros = []
    for col in pvals_adj.columns: 
        if sum(pvals_adj[col]) == 0: 
            zeros.append(col)
    temp = pd.DataFrame(dict(zip(pvals_adj.columns, [1 if val in zeros else 0 for val in pvals_adj.columns])), index = ["unassigned"])
    print(temp.columns)
    print(temp.index)
    # temp = pd.DataFrame(dict(zip(zeros, [1] * len(zeros))), index = ["unassigned"])
    pvals_adj = pd.concat([pvals_adj, temp])
#     pvals_adj.index = pvals_adj.index.set_names('query_cluster', level=0)
#     pvals_adj.columns.name = ["query_cluster", None]
#     pvals_adj.index.name = "ref_cluster"
    print(pvals_adj.columns)
    print(pvals_adj.index)
    return pvals_adj

def reorder_FRmatch(df, axis = 1): 
#     # first try
#     temp = pd.DataFrame(df.head(df.shape[0]-1).sum(0))
#     temp[0] = [1 if val > 0 else 2 for val in temp[0]]
#     temp = temp.sort_values(0, ascending = True)
#     df = df[temp.index]
    
#     # second try
#     todo = list(df.columns)
#     order = []
#     for index, row in df.iterrows(): 
#         for col in todo: 
#             if row[col] != 0: 
#                 order.append(col)
#                 todo.remove(col)
#         if len(todo) == 0: break
#     order.extend(todo)
#     df = df[order]
    
    # third try
    if axis == 0: 
        df = df.sort_values(by = list(df.columns), axis = 0, ascending = False)
    elif axis == 1: 
        df = df.sort_values(by = list(df.index), axis = 1, ascending = False)
    else: print("Error: `axis` must be in [0, 1]")
    
    return df
